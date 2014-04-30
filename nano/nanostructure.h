#ifndef NANOSTRUCTURE_H
#define NANOSTRUCTURE_H

#include "../math/cfunction.h"
#include "../math/function.h"
#include "../dmft/dos.h"
#include "../math/physics.h"
#include "../config/configuration.h"
#include "layerparameters.h"

#include <map>

namespace libconfig {
    class Config;
}

namespace mpi {
    class OpenMPI;
}

namespace nano {
class NanoStructure;
struct IntegrationParametersG {
    NanoStructure* structure;
    int nLayer;
};
struct IntegrationParametersSigmaXX {
    NanoStructure* structure;
    int alpha;
    int beta;
};

struct IntegrationParametersPiXY {
    NanoStructure* structure;
    int alpha;
    int gamma;
    int delta;
};

class NanoStructure
{
public:
    NanoStructure();
    NanoStructure(std::string filename);
    NanoStructure(int n);
    ~NanoStructure();

    void configure(config::Configuration& config);

    void solve();
    void solve(mpi::OpenMPI& mpi);
    std::complex<double> left(int nLayer, int nOmega, double epsk);
    std::complex<double> right(int nLayer, int nOmega, double epsk);

    void precalc_L(int nLayer, int nOmega, double epsk);
    void precalc_R(int nLayer, int nOmega, double epsk);
    void precalc_dL(int nLayer, int nOmega, double epsk);
    void precalc_dR(int nLayer, int nOmega, double epsk);

    double Pi(int alpha, int gamma, int delta, int nOmega, double epsk);

    std::complex<double> G(int nLayer, int nOmega, double epsk);
    std::complex<double> F(int nLayer, int nOmega, double epsk);
    void F_local(int nLayer, math::CFunction& fLocal);
    void G_local(int nLayer, math::CFunction& gLocal);
    void F_local(int nLayer, math::CFunction& fLocal, mpi::OpenMPI& mpi);
    void G_local(int nLayer, math::CFunction& gLocal, mpi::OpenMPI& mpi);
    std::complex<double> GOff(int alpha, int beta, int nOmega, double epsk);

    //ECR
    void updatePotentials();
    void initializePotentials();
    void updatePotentials(mpi::OpenMPI& mpi);
    double calculateOccupationFromG(int layer, double mu = 0.0);
    double totalExcessCharge(double mu);
    double adjustMu();

    // conductivities
    double calculateSigmaXX(int alpha, int beta);
    double calculatePiXY(int alpha, int gamma, int delta);
    double calculateSigmaXX(int alpha, int beta, mpi::OpenMPI& mpi);
    double calculatePiXY(int alpha, int gamma, int delta, mpi::OpenMPI& mpi);

    //consistency questions
    bool areAllSelfEnergySizesEqual();
    bool isScalingCompatibleWithSelfEnergies();

    // convenience
    void setMaxIterations(int maxIterations) { m_maxIterations = maxIterations; }
    void setSelfEnergies(const math::CFunction& selfEnergy);
    void setHoppings(double tL, double tR);
    void setUs(double U);
    bool isConverged();
    inline double getT() const { return m_T; }

    void setAlphaV(double alphaV) { m_alphaV = alphaV; }
    void setCalculateRhoXX(bool calc) { m_calculateRhoXX = calc; }
    void setCalculateRhoXY(bool calc) { m_calculateRhoXY = calc; }
    void setVerbose(bool verbose) { m_verbose = verbose; }
    void setSymmetricNano(bool symmetric) { m_symmetricNano = symmetric; }
    bool isSymmetricNano() { return m_symmetricNano; }

    void setMu(int nLayer, double mu) { m_mu[nLayer] = mu; }
    void setSelfEnergy(int nLayer, const math::CFunction& S) { m_selfEnergies[nLayer] = S; }
    void setU(int nLayer, double U) { m_U[nLayer] = U; }
    void setDoECR(bool ecr) { m_doECR = ecr; }
    //io
    void writeBinary(std::string directory);
    void readBinary(std::string directory);
    void showInfo();
protected:
    double doNRG(int layer);
    void calculateEffectiveMedia(mpi::OpenMPI& mpi);
    void allocate(int N);
    void free();
    typedef std::map<std::string, nano::LayerParameters> tPMap;
    tPMap parseStructure(const libconfig::Config& structure);

    int m_nLayers;
    math::Function* m_n;
    math::CFunction* m_selfEnergies;
    math::CFunction* m_greensFunctions;
    math::CFunction* m_fFunctions;
    math::CFunction* m_effectiveMedia;

    // physical parameters
    // -- on-site U
    double* m_U;
    double* m_mu;
    double m_T;

    // electronic charge reconstruction
    bool m_doECR;
    double* m_rho;
    double* m_rhoBackground;
    double* m_V;
    double m_rhoBar;

    // symmetries
    bool m_symmetricNano;

    // -- hopping
    double* m_tL;
    double* m_tR;


    //numerical parameters
    double* m_similarities;
    double m_delta;
    double m_scaling;
    double m_alphaV;
    double m_tolerance;
    int m_iteration;
    int m_maxIterations;

    //post processing
    bool m_calculateRhoXX;
    bool m_calculateRhoXY;

    bool m_verbose;

    //temp space
    std::complex<double> *m_L;
    std::complex<double> *m_dL;
    std::complex<double> *m_R;
    std::complex<double> *m_dR;
    std::complex<double> *m_temp;

    config::Configuration* m_configuration;


    // wrapper for fermi function
    inline static std::complex<double> fermiWrapper(double omega, std::complex<double> fomega, void* parameters) {
        double * beta = static_cast<double*>(parameters);
        return math::physics::fermi(*beta, omega);
    }

    //integrands
    inline static double gReal(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        const IntegrationParametersG* parameters = static_cast<const IntegrationParametersG*>(params);
         return dmft::dos::rho2(epsk) * parameters->structure->G(parameters->nLayer, nOmega, epsk ).real();
    }
    inline static double gImag(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        const IntegrationParametersG* parameters = static_cast<const IntegrationParametersG*>(params);
        return dmft::dos::rho2(epsk) * parameters->structure->G(parameters->nLayer, nOmega, epsk ).imag();
    }
    inline static double fReal(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        const IntegrationParametersG* parameters = static_cast<const IntegrationParametersG*>(params);
         return dmft::dos::rho2(epsk) * parameters->structure->F(parameters->nLayer, nOmega, epsk ).real();
    }
    inline static double fImag(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        const IntegrationParametersG* parameters = static_cast<const IntegrationParametersG*>(params);
        return dmft::dos::rho2(epsk) * parameters->structure->F(parameters->nLayer, nOmega, epsk ).imag();
    }

    //conductivity
    //    inline static double sigmaXX(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
    //    const IntegrationParametersSigmaXX* parameters = static_cast<const IntegrationParametersSigmaXX*>(params);
    //    return std::pow(parameters->structure->GOff(parameters->alpha, parameters->beta, nOmega, epsk ).imag(),2) * math::physics::dfermi( 1.0 / parameters->structure->getT(), epsk);
    //}

    //conductivity
    inline static double sigmaXX(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
      const IntegrationParametersSigmaXX* parameters = static_cast<const IntegrationParametersSigmaXX*>(params);
      return std::pow(parameters->structure->GOff(parameters->alpha, parameters->beta, nOmega, epsk ).imag(),2) * dmft::dos::rhoXX(epsk);
    }

    inline static double piXY(double epsk, int nOmega, double omega, std::complex<double> fomega, void * params) {
        const IntegrationParametersPiXY* parameters = static_cast<const IntegrationParametersPiXY*>(params);
        return parameters->structure->Pi(parameters->alpha, parameters->gamma, parameters->delta, nOmega, epsk);
    }

};
}

#endif // NANOSTRUCTURE_H
