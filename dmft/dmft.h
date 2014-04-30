#ifndef DMFT_H
#define DMFT_H

#include "../math/cfunction.h"
#include "../utils/complex.h"
#include "dos.h"

namespace config {
    class Configuration;
}

namespace dmft {
class DMFT
{
public:
    DMFT();

    void configure(config::Configuration& configuration);

    //getters and setters
    void setU(double U) { m_U = U; }
    void setMu(double mu) { m_mu = mu; }
    void setT(double T) { m_T = T; }
    void setDelta(double delta) { m_delta = delta; }
    void setTolerance(double tolerance) { m_tolerance = tolerance; }
    void setMaxIterations(int maxIterations) { m_maxIterations = maxIterations; }
    inline double getT() { return m_T; }
    inline double getU() { return m_U; }
    inline double getMu() { return m_mu; }
    inline double getDelta() { return m_delta; }
    inline double getScaling() { return m_scaling; }
    void setStartWithDelta(bool startWithDelta) { m_startWithDelta = startWithDelta; }
    void setScaling(double scaling) { m_scaling = scaling; }
    math::CFunction getSelfEnergy() { return m_selfEnergy; }
    math::CFunction getGreensFunction() { return m_greenFunction; }
    math::CFunction getHybridizationFunction() { return m_effectiveMedium; }
    void setSelfEnergy(const math::CFunction& selfEnergy) { m_selfEnergy = selfEnergy; }
    void setGreensFunction(const math::CFunction& G) { m_greenFunction = G; }
    void setHybridizationFunction(const math::CFunction& delta) { m_effectiveMedium = delta; }
    void solve(std::string outputDirectory = "");
    ~DMFT();
protected:

    inline static double gReal(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        DMFT* dmft = static_cast<DMFT*>(params);
        std::complex<double> G = 1.0 / (omega + I * dmft->getDelta() - fomega - (epsk - dmft->getU()/2.0 - dmft->getMu()));
        return G.real() * dos::rho3(epsk);
    }

    inline static double gImag(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        DMFT* dmft = static_cast<DMFT*>(params);
        std::complex<double> G = 1.0 / (omega + I * dmft->getDelta() - fomega - (epsk - dmft->getU()/2.0 - dmft->getMu()));
        return G.imag() * dos::rho3(epsk);
    }

    inline static double fReal(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        DMFT* dmft = static_cast<DMFT*>(params);
        std::complex<double> G = 1.0 / (omega + I * dmft->getDelta() - fomega - (epsk - dmft->getU()/2.0 - dmft->getMu()));
        std::complex<double> F = fomega * G + 1.0;
        return F.real() * dos::rho3(epsk);
    }

    inline static double fImag(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        DMFT* dmft = static_cast<DMFT*>(params);
        std::complex<double> G = 1.0 / (omega + I * dmft->getDelta() - fomega - (epsk - dmft->getU()/2.0 - dmft->getMu()));
        std::complex<double> F = fomega * G + 1.0;
        return F.imag() * dos::rho3(epsk);
    }

    void calculateGFCF();
    void calculateGFInt();
    void calculateFFInt();

    //helper functions
    double Im_sqr(double Re, double Im);
    double Re_sqr(double Re, double Im);
    double Im_sqrt(double Re, double Im);
    double Re_sqrt(double Re, double Im);
    double Im_calc(double Re, double Im);
    double Re_calc(double Re, double Im);

    double m_U;
    double m_T;
    double m_mu;
    double m_delta;

    double m_scaling;
    double m_tolerance;
    int m_maxIterations;

    bool m_startWithDelta;

    math::CFunction m_selfEnergy;
    math::CFunction m_greenFunction;
    math::CFunction m_fFunction;
    math::CFunction m_effectiveMedium;

    config::Configuration* m_configuration;
};
}
#endif // DMFT_H
