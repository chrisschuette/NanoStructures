#ifndef DMFTAB_H
#define DMFTAB_H

#include "../math/cfunction.h"
#include "dos.h"

namespace config {
    class Configuration;
}

namespace dmft {
class DMFTAB
{
public:
    DMFTAB();
    ~DMFTAB();

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
    math::CFunction getSelfEnergyA() const { return m_selfEnergyA; }
    math::CFunction getGreensFunctionA() const { return m_greenFunctionA; }
    math::CFunction getHybridizationFunctionA() const { return m_effectiveMediumA; }
    math::CFunction getSelfEnergyB() const { return m_selfEnergyB; }
    math::CFunction getGreensFunctionB() const { return m_greenFunctionB; }
    math::CFunction getHybridizationFunctionB() const { return m_effectiveMediumB; }
    double getOccupation() const { return m_n; }
    double getMagnetization() const { return m_m; }
    void setSelfEnergyA(const math::CFunction& selfEnergy) { m_selfEnergyA = selfEnergy; }
    void setGreensFunctionA(const math::CFunction& G) { m_greenFunctionA = G; }
    void setHybridizationFunctionA(const math::CFunction& delta) { m_effectiveMediumA = delta; }
    void setSelfEnergyB(const math::CFunction& selfEnergy) { m_selfEnergyB = selfEnergy; }
    void setGreensFunctionB(const math::CFunction& G) { m_greenFunctionB = G; }
    void setHybridizationFunctionB(const math::CFunction& delta) { m_effectiveMediumB = delta; }

    void solve(std::string outputDirectory = "");
protected:
    inline static double gAReal(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        std::complex<double> z(omega, dmft->getDelta());
        std::complex<double> XAup = z + dmft->getMu() + dmft->getU()/2.0 - fomegaA;
        std::complex<double> XBup = z + dmft->getMu() + dmft->getU()/2.0 - fomegaB;
        return (XBup / (XAup * XBup - epsk * epsk)).real() * dos::rho3(epsk);
    }

    inline static double gAImag(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        std::complex<double> z(omega, dmft->getDelta());
        std::complex<double> XAup = z + dmft->getMu() + dmft->getU()/2.0 - fomegaA;
        std::complex<double> XBup = z + dmft->getMu() + dmft->getU()/2.0 - fomegaB;
        return (XBup / (XAup * XBup - epsk * epsk)).imag() * dos::rho3(epsk);
    }

    inline static double gBReal(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        std::complex<double> z(omega, dmft->getDelta());
        std::complex<double> XAdown = z + dmft->getMu() + dmft->getU()/2.0 - fomegaB;
        std::complex<double> XBdown = z + dmft->getMu() + dmft->getU()/2.0 - fomegaA;
        return (XBdown / (XAdown * XBdown - epsk * epsk)).real() * dos::rho3(epsk);
    }

    inline static double gBImag(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        std::complex<double> z(omega, dmft->getDelta());
        std::complex<double> XAdown = z + dmft->getMu() + dmft->getU()/2.0 - fomegaB;
        std::complex<double> XBdown = z + dmft->getMu() + dmft->getU()/2.0 - fomegaA;
        return (XBdown / (XAdown * XBdown - epsk * epsk)).imag() * dos::rho3(epsk);
    }

    inline static double fAReal(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        return 0;
    }

    inline static double fAImag(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        return 0;

    }

    inline static double fBReal(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        return 0;

    }

    inline static double fBImag(double epsk, int nOmega, double omega, std::complex<double> fomegaA, std::complex<double> fomegaB, void * params) {
        DMFTAB* dmft = static_cast<DMFTAB*>(params);
        return 0;

    }


    void calculateGFAInt();
    void calculateGFBInt();
    void calculateFFAInt();
    void calculateFFBInt();



    double m_U;
    double m_T;
    double m_mu;
    double m_delta;

    double m_scaling;
    double m_tolerance;
    int m_maxIterations;

    bool m_startWithDelta;

    math::CFunction m_selfEnergyA;
    math::CFunction m_selfEnergyB;
    math::CFunction m_greenFunctionA;
    math::CFunction m_greenFunctionB;
    math::CFunction m_fFunctionA;
    math::CFunction m_fFunctionB;
    math::CFunction m_effectiveMediumA;
    math::CFunction m_effectiveMediumB;

    double m_n;
    double m_m;

    config::Configuration* m_configuration;
};
}

#endif // DMFTAB_H
