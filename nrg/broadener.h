#ifndef BROADENER_H
#define BROADENER_H

#include "../math/cfunction.h"
#include "../config/configuration.h"

namespace nrg {
class Broadener
{
public:
    Broadener();
    Broadener(const Broadener& orig);
    Broadener& operator=(const Broadener& orig);
    virtual ~Broadener();

    virtual void configure(config::Configuration& configuration);

    void addExcitation(double omega, double weight);
    void init();

    virtual math::CFunction broaden() { return math::CFunction(); }
    virtual Broadener* clone() { return new Broadener(*this); }
    virtual void showInfo();
    virtual void showStatistics();

    // getters and setters
    void setPolesPerDecade(int polesPerDecade) { m_polesPerDecade = polesPerDecade; }
    void setPeakWidth(double peakWidth) { m_peakWidth = peakWidth; }
    void setTemperature(double temperature) { m_temperature = temperature; }

protected:
    int mapFrequencyToGrid(double frequency);

    double m_gamma;
    double m_temperature;
    double m_lambda;
    int m_polesPerDecade;
    double m_weightInSpectrum;
    double m_w0;
    double m_w0Weight;
    int m_w0PeakCount;
    int m_peakCount;
    double m_omegaMin;
    double m_omegaMax;
    double m_zeroWeight;
    double m_logLambda;
    int m_lMax;
    double * m_weightsNegative;
    double * m_weightsPositive;
    double m_negativeWeight;
    double m_positiveWeight;
    double m_peakWidth;
};
}
#endif // BROADENER_H
