#ifndef NRG_H
#define NRG_H

#include "../config/configuration.h"
#include "chainprovider.h"
#include "hilbertspacetable.h"
#include "broadener.h"

namespace nrg {
class NRG
{
    enum STATES {KEPT, DISCARDED, BOTH};

public:
    NRG(ChainProvider& chainProvider, Broadener& broadener);
    void configure(config::Configuration& configuration);

    // setters and getters -- physical parameters
    double getU() const { return m_U; }
    void setU(double U) { m_U = U; }

    double getEpsF() const { return m_epsF; }
    void setEpsF(double epsF) { m_epsF = epsF; }

    double getTemperature() const { return m_temperature; }
    void setTemperature( double T ) { m_temperature = T; }

    //setters and getters -- numerical parameters
    double getClusterEnergy() const { return m_clusterEnergy; }
    void setClusterEnergy(double clusterEnergy) { m_clusterEnergy = clusterEnergy; }

    double getEnergyCutOff() const { return m_energyCutOff; }
    void setEnergyCutOff(double energyCutOff) { m_energyCutOff = energyCutOff; }

    int getMaxHilbertSpaceDimension() { return m_maxHSdimension; }
    void setMaxHilbertSpaceDimension(int maxHSdimension) { m_maxHSdimension = maxHSdimension; }

    int getMaxIterations() { return m_maxIterations; }
    void setMaxIterations(int maxIterations) { m_maxIterations = maxIterations; }

    inline double getOccupationUp() { return m_n_Up; }
    inline double getOccupationDown() { return m_n_Down; }
    inline double getOccupation() { return m_n_Up + m_n_Down; }
    inline double getMagnetization() { return m_n_Up - m_n_Down; }

    // nrg
    void init();
    // - determine m_maxIterations from m_temperature (if necessary)
    // - call build chain
    // - call setupInitialState();
    void setupInitialState();
    void setupHamiltonian(int iteration); // iteration >= 0
    void truncateStates(int iteration);

    void propagateChainOperatorElementsUp(int iteration);
    void propagateChainOperatorElementsDown(int iteration);

    void propagateLocalMatrixElementUp(int iteration);
    void propagateLocalMatrixElementDown(int iteration);
    void propagateLocalMatrixElementUp2(int iteration);
    void propagateLocalMatrixElementDown2(int iteration);
    void createPolesG_Up(int iteration);
    void createPolesG_Down(int iteration);
    void createPolesF_Up(int iteration);
    void createPolesF_Down(int iteration);

    void solve(bool silent = false);
    void solve_symmetric_SZ(bool silent = false);

    void builDM();

    // show parameters
    void showInfo();

    // clean-up code
    void deleteChainOperatorElements(int iteration);
    void deleteTransformationMatrices(int iteration);
    void deleteDensityMatrices(int iteration);
    void deleteImpurityMatrixElements(int iteration);
    ~NRG();

    void getSelfEnergy(math::CFunction& SUp, math::CFunction& SDown) { SUp = m_F_S_Up; SDown = m_F_S_Down; }
    void getSelfEnergy(math::CFunction& S) { S = m_F_S_Up; }

    void getGreensFunction(math::CFunction& GUp, math::CFunction& GDown) { GUp = m_F_G_Up; GDown = m_F_G_Down; }
    void getGreensFunction(math::CFunction& G) { G = m_F_G_Up; }

    void getFFunction(math::CFunction& FUp, math::CFunction& FDown) { FUp = m_F_F_Up; FDown = m_F_F_Down; }
    void getFFunction(math::CFunction& F) { F = m_F_F_Up; }

protected:
    ChainProvider& m_chainProvider;
    Broadener& m_broadener;
    Broadener* m_G_Up;
    Broadener* m_G_Down;
    Broadener* m_F_Up;
    Broadener* m_F_Down;

    math::CFunction m_F_S_Up;
    math::CFunction m_F_S_Down;
    math::CFunction m_F_G_Up;
    math::CFunction m_F_G_Down;
    math::CFunction m_F_F_Up;
    math::CFunction m_F_F_Down;

    double m_n_Up;
    double m_n_Down;

    double m_epsF;
    double m_U;
    double m_temperature;
    int m_maxIterations;
    std::vector<double> m_energies;

    HilbertSpaceTable m_hilbertSpaces;

    // numerical parameters
    double m_clusterEnergy;
    double m_energyCutOff;
    int m_maxHSdimension;
    int m_nFirstTruncated;

    // helper fields
    int signLME[4];
    int dQ[4];
    int dSz[4];

};
}
#endif // NRG_H
