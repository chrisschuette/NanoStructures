#include "bullabroadener.h"
#include <iostream>

namespace nrg {
namespace broadening {


double Lorentzian(double freq, double center_freq, double width)
{
    return M_PI*width/( M_PI * ( freq - center_freq ) * ( freq - center_freq ) + width * width );
}

double LogGaussian(double freq, double center_freq, double width)
{
    return M_PI*exp(-width*width/4.0)/( width * fabs(center_freq) * sqrt(M_PI) ) * exp( -std::pow(log(fabs(freq))-log(fabs(center_freq)),2.0)/(width * width)  );
}

inline double f(double x, double b)
{
    return exp(-b*b/4.0)/(b*sqrt(M_PI))*exp(-(log(x)/b)*(log(x)/b));
}

inline double KKTLogGaussian_integrand(double x, double z, double b,double delta)
{
    return f(z/(1.0-z),b)*(x-z/(1.0-z))/( std::pow(delta * (z-1.0),2.0)+std::pow((x-z/(1.0-z))*(z-1.0),2.0 ) );
}

inline double KKTLorentzian(double x, double x0, double b,double delta)
{
    return (x0-x)/((std::pow((b-delta),2.0) + std::pow(x-x0,2.0) )   );
}

double KKTLogGaussian(double x,double width)
{

    double integral = 0.0;

    int steps = 1000;

    double xR = 0.0;
    double xL = 0.0;

    double integrandR = KKTLogGaussian_integrand(x,0.0,width,0.01);
    double integrandL = integrandR;
    for(int i = 0; i < steps-1; i++)
    {
        xL = xR;
        xR = (i+1.0)/1000.0;
        integrandL = integrandR;
        integrandR = KKTLogGaussian_integrand(x,xR,width,0.01);
        integral += (xR - xL)*(integrandL + integrandR)/2.0;
    }
    return integral;
}

BullaBroadener::BullaBroadener()
    : Broadener()
{
}

BullaBroadener::BullaBroadener(const BullaBroadener& orig)
    : Broadener(orig) {

}

BullaBroadener& BullaBroadener::operator=(const BullaBroadener& orig) {
    Broadener::operator =(orig);
    return *this;
}

void BullaBroadener::showInfo() {

}

math::CFunction BullaBroadener::broaden() {
    m_cachePositive = math::Function(2 * m_lMax - 1);
    for(int iminusj = -(m_lMax-1); iminusj <= (m_lMax-1); iminusj++)
    {
            double w = std::pow((double) m_lambda, (double) -iminusj);
            double value = KKTLogGaussian(w,m_peakWidth);
            m_cachePositive.set(iminusj + (m_lMax-1),w,value);
    }
    m_cacheNegative = math::Function(2 * m_lMax - 1);
    for(int iminusj = -(m_lMax-1); iminusj <= (m_lMax-1); iminusj++)
    {
            double w = -std::pow((double) m_lambda, (double) -iminusj);
            double value = KKTLogGaussian(w,m_peakWidth);
            m_cacheNegative.set(iminusj + (m_lMax-1),w,value);
    }
    double b_width = 0.7;

double b_LogGauss = b_width;
double b_Lorentzian = 2 * b_width * m_temperature;

int refinement = 10;

std::cout << "Grid: Lambda: " << m_lambda << " L_MAX: " << m_lMax << std::endl;



std::cout << "Temperature: " << m_temperature << std::endl;

double w0 = 4*m_temperature;

    std::cout << "Broadening" << std::endl;
    long int nmax = m_lMax;
    math::CFunction broadened(2 * nmax);
    for(long int j = 0; j < nmax; j++)
    {
            double freq = -std::pow((double) m_lambda, (double) -j);
            double weight = 0.0;
            double real_weight = 0.0;
    //	cout << freq << nendl;


                    //negative freq

                    //negative peaks
                    for(long int i = 0; i <= m_lMax; i++)
                    {
                            double peak_freq = -m_omegaMax * std::pow((double) m_lambda, (double) -i);
                            double peak_weight = m_weightsNegative[i];
                            if(std::abs(peak_freq) < w0)
                            {
                                    //lorentzian
                                    weight += peak_weight * Lorentzian(freq,peak_freq,b_Lorentzian);
                                    real_weight += peak_weight * KKTLorentzian(freq,peak_freq,b_Lorentzian,0.01);
                            }
                            else
                            {
                                    //log gaussian
                                    weight += peak_weight * LogGaussian(freq,peak_freq,b_LogGauss);
                                    real_weight += peak_weight * m_cachePositive.getValue(j-refinement*i + m_lMax - 1)/peak_freq;
                            }
                    }

                    //zero peak
                    //weight  += b_Lorentzian * Lorentzian(freq,0.0,b_Lorentzian);



                    //positive peaks
                    for(long int i = 0; i <= m_lMax; i++)
                    {
                            double peak_freq = m_omegaMax * std::pow((double) m_lambda, (double) -i);
                            double peak_weight = m_weightsPositive[i];
                            if(std::abs(peak_freq) < w0)
                            {
                                    //lorentzian
                                    weight += peak_weight * Lorentzian(freq,peak_freq,b_Lorentzian);
                                    real_weight += peak_weight * KKTLorentzian(freq,peak_freq,b_Lorentzian,0.01);
                            }
                            else
                            {
                                    //log gaussian
                                    real_weight += peak_weight * m_cacheNegative.getValue(j-refinement*i + m_lMax - 1)/peak_freq;
                            }
                    }

                    broadened.set(j,freq, std::complex<double>(real_weight, -weight));
    }

    for(long int j = nmax-1; j >= 0; j--)
    {
            double freq = std::pow((double) m_lambda, (double) -j);
            double weight = 0.0;
            double real_weight = 0.0;
    //	cout << freq << nendl;

                    //positive freq


                    //negative peaks
                    for(long int i = 0; i <= m_lMax; i++)
                    {
                            double peak_freq = -m_omegaMax * std::pow((double) m_lambda, (double) -i);
                            double peak_weight = m_weightsNegative[i];
                            if(std::abs(peak_freq) < w0)
                            {
                                    //lorentzian
                                    weight += peak_weight * Lorentzian(freq,peak_freq,b_Lorentzian);
                                    real_weight += peak_weight * KKTLorentzian(freq,peak_freq,b_Lorentzian,0.01);
                            }
                            else
                            {
                                    //log gaussian
                                    real_weight += peak_weight * m_cacheNegative.getValue(j-refinement*i + m_lMax - 1)/peak_freq;
                            }
                    }

                    //zero peak
                    //weight  += b_Lorentzian * Lorentzian(freq,0.0,b_Lorentzian);



                    //positive peaks
                    for(long int i = 0; i <= m_lMax; i++)
                    {
                            double peak_freq = m_omegaMax * std::pow((double) m_lambda, (double) -i);
                            double peak_weight = m_weightsPositive[i];
                            if(std::abs(peak_freq) < w0)
                            {
                                    //lorentzian
                                    weight += peak_weight * Lorentzian(freq,peak_freq,b_Lorentzian);
                                    real_weight += peak_weight * KKTLorentzian(freq,peak_freq,b_Lorentzian,0.01);
                            }
                            else
                            {
                                    //log gaussian
                                    weight +=  peak_weight * LogGaussian(freq,peak_freq,b_LogGauss);
                                    real_weight += peak_weight * m_cachePositive.getValue(j-refinement*i + m_lMax - 1)/peak_freq;
                            }
                    }



    broadened.set(nmax + nmax-1-j,freq, std::complex<double>(real_weight, -weight));
    }

return broadened;
}

}
}
