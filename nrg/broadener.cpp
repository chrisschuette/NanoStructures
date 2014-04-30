#include "broadener.h"
#include <cmath>
#include <iostream>
#include <string.h>

namespace nrg {
Broadener::Broadener()
    : m_temperature(1e-20)
    , m_lambda(0)
    , m_polesPerDecade(50)
    , m_weightInSpectrum(0)
    , m_w0(0)
    , m_w0Weight(0)
    , m_w0PeakCount(0)
    , m_peakCount(0)
    , m_omegaMin(1e-20)
    , m_omegaMax(1.0)
    , m_zeroWeight(0)
    , m_logLambda(0)
    , m_lMax(0)
    , m_weightsNegative(0)
    , m_weightsPositive(0)
    , m_negativeWeight(0)
    , m_positiveWeight(0)
    , m_peakWidth(0.7)
{
}

Broadener::Broadener(const Broadener& orig)
    : m_temperature(orig.m_temperature)
    , m_lambda(orig.m_lambda)
    , m_polesPerDecade(orig.m_polesPerDecade)
    , m_weightInSpectrum(orig.m_weightInSpectrum)
    , m_w0(orig.m_w0)
    , m_w0Weight(orig.m_w0Weight)
    , m_w0PeakCount(orig.m_w0PeakCount)
    , m_peakCount(orig.m_peakCount)
    , m_omegaMin(orig.m_omegaMin)
    , m_omegaMax(orig.m_omegaMax)
    , m_zeroWeight(orig.m_zeroWeight)
    , m_logLambda(orig.m_logLambda)
    , m_lMax(orig.m_lMax)
    , m_negativeWeight(orig.m_negativeWeight)
    , m_positiveWeight(orig.m_positiveWeight)
    , m_peakWidth(orig.m_peakWidth)
{
    if(orig.m_weightsNegative != 0) {
        m_weightsNegative = new double[m_lMax+1];
        memcpy(m_weightsNegative, orig.m_weightsNegative, sizeof(double) * (m_lMax+1));
    }
    else
        m_weightsNegative = 0;
    if(orig.m_weightsPositive != 0) {
        m_weightsPositive = new double[m_lMax+1];
        memcpy(m_weightsPositive, orig.m_weightsPositive, sizeof(double) * (m_lMax+1));
    }
    else
        m_weightsPositive = 0;
}

Broadener& Broadener::operator=(const Broadener& orig) {
    m_temperature = orig.m_temperature;
    m_lambda = orig.m_lambda;
    m_polesPerDecade = orig.m_polesPerDecade;
    m_weightInSpectrum = orig.m_weightInSpectrum;
    m_w0 = orig.m_w0;
    m_w0Weight = orig.m_w0Weight;
    m_w0PeakCount = orig.m_w0PeakCount;
    m_peakCount = orig.m_peakCount;
    m_omegaMin = orig.m_omegaMin;
    m_omegaMax = orig.m_omegaMax;
    m_zeroWeight = orig.m_zeroWeight;
    m_logLambda = orig.m_logLambda;
    m_lMax = orig.m_lMax;
    m_negativeWeight = orig.m_negativeWeight;
    m_positiveWeight = orig.m_positiveWeight;
    m_peakWidth = orig.m_peakWidth;

    if(m_weightsPositive != 0)
        delete [] m_weightsPositive;
    if(m_weightsNegative != 0)
        delete [] m_weightsNegative;

    if(orig.m_weightsNegative != 0) {
        m_weightsNegative = new double[m_lMax+1];
        memcpy(m_weightsNegative, orig.m_weightsNegative, sizeof(double) * (m_lMax+1));
    }
    else
        m_weightsNegative = 0;
    if(orig.m_weightsPositive != 0) {
        m_weightsPositive = new double[m_lMax+1];
        memcpy(m_weightsPositive, orig.m_weightsPositive, sizeof(double) * (m_lMax+1));
    }
    else
        m_weightsPositive = 0;

    return *this;
}

Broadener::~Broadener() {
    if(m_weightsNegative != 0)
        delete [] m_weightsNegative;
    if(m_weightsPositive != 0)
        delete [] m_weightsPositive;
}

void Broadener::configure(config::Configuration& configuration) {
    try {
        m_polesPerDecade = configuration.getInteger("Broadener.PPD");
    } catch( libconfig::SettingNotFoundException ) {}
    try {
        m_peakWidth = configuration.getDouble("Broadener.peakWidth");
    } catch( libconfig::SettingNotFoundException ) {}
}

void Broadener::init() {
    m_w0 = m_temperature / 2.0;

    m_lambda = exp(log(10.0) / (double) m_polesPerDecade);
    m_logLambda = log(m_lambda);
  //  std::cout << "Discretization parameter for broadening grid: " << m_lambda << std::endl;

    //setup grid
    m_omegaMin = 1e-5 * m_temperature;
    m_omegaMax = 1.0;
    m_lMax = (int) ceil(log(m_omegaMax / m_omegaMin) / log(m_lambda));

  //  std::cout << "Number of Bins in Peak Grid = " << 2 * (m_lMax + 1) << std::endl;

    m_weightsNegative = new double[m_lMax + 1];
    m_weightsPositive = new double[m_lMax + 1];

    //init fields
    m_zeroWeight = 0.0;
    for (long int i = 0; i <= m_lMax; i++) {
        m_weightsNegative[i] = 0.0;
        m_weightsPositive[i] = 0.0;
    }

    m_peakCount = 0;
    m_positiveWeight = 0.0;
    m_negativeWeight = 0.0;
    m_weightInSpectrum = 0.0;
}

void Broadener::addExcitation(double omega, double weight) {
    double absOmega = std::abs(omega);

    m_weightInSpectrum += weight;
    if(absOmega < m_w0)
    {
        m_w0Weight += weight;
        m_w0PeakCount++;
    }
    m_peakCount++;
    if (absOmega < m_omegaMin)
        m_zeroWeight += weight;
    else {
        //select the grid
        double * grid;
        if (omega < 0.0) {
            grid = m_weightsNegative;
            m_negativeWeight += weight;
        } else {
            grid = m_weightsPositive;
            m_positiveWeight += weight;
        }
        //calculate index
        int index = (int) round(log(m_omegaMax / absOmega) / m_logLambda);
        if (index > m_lMax) {
            grid[m_lMax] += weight;
            //ERROR("Peak outside the interval.")
        } else if (index < 0) {
            grid[0] += weight;
        } else {
            if ((index >= 0) && index < (m_lMax + 1))
                grid[index] += weight;
        }
    }

}

int Broadener::mapFrequencyToGrid(double frequency) {
        return (int) round(log(m_omegaMax / fabs(frequency)) / log(m_lambda));
}
void Broadener::showStatistics() {
    std::cout << "============================" << std::endl;
    std::cout << "** Broadener statistics : **" << std::endl;
    std::cout << "============================" << std::endl;
    std::cout << std::endl;
    std::cout << " * logarithmic binning parameters = " << m_lambda << std::endl;
    std::cout << " * minimum frequency = " << m_omegaMin << std::endl;
    std::cout << " * maximum frequency = " << m_omegaMax << std::endl;
    std::cout << " * weight in Spectrum = " << m_weightInSpectrum << std::endl;
    std::cout << " * positive weight = " << m_positiveWeight << std::endl;
    std::cout << " * negative weight = " << m_negativeWeight << std::endl;
    std::cout << std::endl;

}

void Broadener::showInfo() {
    std::cout << "================" << std::endl;
    std::cout << "** Broadener: **" << std::endl;
    std::cout << "================" << std::endl;
    std::cout << std::endl;

    std::cout << " - physical parameters:" << std::endl;
    std::cout << " * temperature = " << m_temperature << std::endl;
    std::cout << " * peakWidth = " << m_peakWidth << std::endl;
    std::cout << " * PPD = " << m_polesPerDecade << std::endl;
    std::cout << std::endl;

}
}
