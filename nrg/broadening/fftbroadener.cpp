#include "fftbroadener.h"
#include "../../math/fourier.h"
#include "../../math/function.h"
#include "../../math/hilberttransformer.h"

#include <iostream>

namespace nrg {
namespace broadening {
FFTBroadener::FFTBroadener()
    :Broadener()
    , m_sqrtPi(sqrt(M_PI))
    , m_epsPadding(1e-3)
{
}

FFTBroadener::FFTBroadener(const FFTBroadener& orig)
    : Broadener(orig), m_sqrtPi(orig.m_sqrtPi), m_epsPadding(orig.m_epsPadding)
{

}

FFTBroadener& FFTBroadener::operator=(const FFTBroadener& orig) {
    Broadener::operator =(orig);
    return *this;
}

FFTBroadener::~FFTBroadener() {

}

int FFTBroadener::getPadding(double eps) {
    return (int) ceil((std::pow(m_peakWidth,2.0)*log(m_lambda) -
      2*m_peakWidth*m_gamma*log(m_lambda) +
      m_peakWidth*sqrt(std::pow(m_peakWidth,2) - 4.0 * m_peakWidth*m_gamma -
         4*log(m_peakWidth*sqrt(M_PI)*eps))*log(m_lambda))/
    (2.*std::pow(log(m_lambda),2.0)));
}

math::CFunction FFTBroadener::broaden() {

    std::complex<double> complexZero(0.0,0.0);

    //Ralf
    //m_gamma = m_peakWidth / 2.0;

    //van Delft
    m_gamma = m_peakWidth / 4.0;

    int padding = getPadding(m_epsPadding);
    int N = (m_lMax + 1) + padding;

    // setup frequency grid
    m_frequencies.resize(2 * (m_lMax + 1) );

    //negative frequencies
    for (int i = 0; i < (m_lMax + 1); i++)
        m_frequencies[i] = -std::pow(m_lambda, (double) -i);

    for (int i = 0; i < (m_lMax + 1); i++)
        m_frequencies[m_lMax + 1 + i] = std::pow(m_lambda, (double) -(m_lMax - i));



    std::vector< std::complex<double> > positivePoles;
    positivePoles.resize(N, complexZero);
    for( int i = 0; i < (m_lMax + 1); i++ )
        positivePoles.at(i) = std::complex<double>(m_weightsPositive[i] / m_frequencies[2 * (m_lMax + 1) - 1 - i], 0.0);
    math::fourier::fourierTransform(positivePoles);

    std::vector< std::complex<double> > negativePoles;
    negativePoles.resize(N, complexZero);
    for( int i = 0; i < (m_lMax + 1); i++ )
        negativePoles.at(i) = std::complex<double>(m_weightsNegative[i] / (-m_frequencies[i]), 0.0);
    math::fourier::fourierTransform(negativePoles);

    std::vector< std::complex<double> > response;
    response.resize(N, complexZero);
    for( int i = 0; i < (padding + 1); i++ )
        response.at(i) = std::complex<double>(logGaussian(i), 0.0);
    for( int i = padding + (N - (2 * padding + 1)) + 1; i < N; i++ )
        response.at(i) = std::complex<double>(logGaussian(-(N - i)), 0.0);
    math::fourier::fourierTransform(response);

    std::vector< std::complex<double> > positiveBroadened;
    positiveBroadened.resize(N, complexZero);
    for( int i = 0; i < N; i++ )
        positiveBroadened.at(i) = std::complex<double>((positivePoles.at(i).real() * response.at(i).real() - positivePoles.at(i).imag() * response.at(i).imag())
                , (positivePoles.at(i).real() *  response.at(i).imag() + positivePoles.at(i).imag() * response.at(i).real()) );
    math::fourier::inverseFourierTransform(positiveBroadened);

    std::vector< std::complex<double> > negativeBroadened;
    negativeBroadened.resize(N, complexZero);
    for( int i = 0; i < N; i++ )
        negativeBroadened.at(i) = std::complex<double>((negativePoles.at(i).real() * response.at(i).real() - negativePoles.at(i).imag() * response.at(i).imag())
                , (negativePoles.at(i).real() *  response.at(i).imag() + negativePoles.at(i).imag() * response.at(i).real())  );
    math::fourier::inverseFourierTransform(negativeBroadened);

    math::CFunction temp(m_lMax + 1 + m_lMax + 1);

    //move data into function
    for (int i = 0; i < (m_lMax + 1); i++)
        temp.set(i,m_frequencies[i], std::complex<double>(0.0, negativeBroadened.at(i).real()));

    for (int i = m_lMax; i >= 0; i--)
        temp.set((m_lMax + 1) + m_lMax - i,m_frequencies[m_lMax + 1 + (m_lMax - i)], std::complex<double>(0.0, positiveBroadened.at(i).real()));

    math::CFunction broadened((2 * (m_lMax + 1) ));

    //now what is missing is the peaks below w0

    //go over all frequencies
    for (int i = 0; i < (2 * (m_lMax + 1) ); i++) {
        double omega = m_frequencies[i];
        double con = 0.0;
        if(fabs(m_frequencies[i]) < m_w0) {
        for (int j = 0; j < (2 * (m_lMax + 1) ); j++) {
                double peakPosition = m_frequencies[j];
                    double peakWeight;
                    if (j < (m_lMax + 1)) {
                        peakWeight = m_weightsNegative[j];
                        con += contribution(omega, peakPosition, peakWeight);
                    }  else {
                        peakWeight = m_weightsPositive[2 * (m_lMax + 1) - 1 - j];
                        con += contribution(omega, peakPosition, peakWeight);
                    }
            }
        }
        broadened.set(i,m_frequencies[i],  M_PI * (std::complex<double>(0.0, 1.0) * con * (1.0 - getCrossOver(m_frequencies[i]))  + getCrossOver(m_frequencies[i]) * temp.getValue(i) )  );

    }


    math::Function posPoles(m_lMax+1);
    math::Function negPoles(m_lMax+1);

    for(int i = m_lMax; i >= 0; i--)
        posPoles.set(m_lMax-i,m_frequencies[i], broadened.getValue(i).imag() );

    for(int i = m_lMax+1; i < (2 * (m_lMax + 1) ); i++)
        negPoles.set(i-(m_lMax+1),m_frequencies[i], broadened.getValue(i).imag() );

    math::HilbertTransformer ht;
    ht.setPoles(posPoles,negPoles);
    ht.init();
    ht.broaden();


    for(unsigned int i = 0; i < ht.posHT.size() ; i++)
    {
        std::complex<double> temp;
        temp = broadened.getValue(m_lMax-i);
        temp.real() = ht.posHT[i];
        broadened.setValue(m_lMax-i,temp);
    }

    for(unsigned int i = 0; i < ht.negHT.size() ; i++)
    {
        std::complex<double> temp;
        temp = broadened.getValue(m_lMax+i+1);
        temp.real() = ht.negHT[i];
         broadened.setValue(m_lMax+i+1,temp);
    }

    int counter = 0;
    for( int i = 0; i < broadened.getSize(); i++)
        if(fabs(broadened.getArgument(i)) > (m_w0 / 100) )
            counter++;


    //TODO: maybe we should throw away the data here which is not trustworthy...
    math::CFunction broadened2(counter);

    counter = 0;
    for( int i = 0; i < broadened.getSize(); i++)
        if(fabs(broadened.getArgument(i)) > (m_w0 / 100) ) {
            broadened2.set(counter,broadened.getArgument(i), broadened.getValue(i));
            counter++;
        }
    return broadened2;
}

double FFTBroadener::getCrossOver(double frequency) {
    double arg = fabs(frequency / m_w0);
    if (arg >= 1.0)
        return 1.0;
    else
        return exp(-(std::pow(log(arg) / (m_peakWidth), 2.0)));
}

double FFTBroadener::gaussian(double frequency) {
    double b = m_w0;
    double a = (frequency / b);
    return 1.0 / (m_sqrtPi * b) * exp(-a*a);
}

double FFTBroadener::contribution(double omega, double peakPosition, double peakWeight) {
    return gaussian(omega - peakPosition) * peakWeight;// * (1.0 - getCrossOver(peakPosition));
}

double FFTBroadener::logGaussian(int j) {
    return exp(-(m_peakWidth*(-m_peakWidth/4. + m_gamma)) -
     std::pow(-m_peakWidth/2. + m_gamma + (j*log(m_lambda))/m_peakWidth,2.0))/
   (m_peakWidth*sqrt(M_PI));

}

void FFTBroadener::showInfo() {
    Broadener::showInfo();
}

}
}
