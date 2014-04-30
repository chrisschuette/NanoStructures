#include "dmftab.h"
#include "../config/configuration.h"
#include "../nrg/chain/hybridizationprovider.h"
#include "../nrg/broadening/fftbroadener.h"
#include "../nrg/nrg.h"

#include "../math/integralAB.h"
#include "../utils/utils.h"
#include "../utils/complex.h"

#include <iostream>

namespace dmft {
DMFTAB::DMFTAB()
    : m_U(0.0)
    , m_T(1e-8)
    , m_mu(0.0)
    , m_delta(1e-3)
    , m_scaling(20.0)
    , m_tolerance(0.9995)
    , m_maxIterations(100)
    , m_startWithDelta(false)
    , m_configuration(0)
{
}

DMFTAB::~DMFTAB() {

}


void DMFTAB::configure(config::Configuration& configuration) {
    m_configuration = &configuration;
    try {
        m_U = configuration.getDouble("DMFTAB.U");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_mu = configuration.getDouble("DMFTAB.mu");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_T = configuration.getDouble("DMFTAB.temperature");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_delta = configuration.getDouble("DMFTAB.delta");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_scaling = configuration.getDouble("DMFTAB.scaling");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_tolerance = configuration.getDouble("DMFTAB.tolerance");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_maxIterations = configuration.getInteger("DMFTAB.maxIterations");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_startWithDelta = configuration.getBool("DMFTAB.startWithDelta");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        std::string initialDA = configuration.getString("DMFTAB.initialDA");
        std::string initialDB = configuration.getString("DMFTAB.initialDB");
        math::CFunction DA;
        math::CFunction DB;
        DA.readBinary(initialDA);
        setHybridizationFunctionA(DA);
        DB.readBinary(initialDB);
        setHybridizationFunctionB(DB);
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        std::string initialSA = configuration.getString("DMFTAB.initialSA");
        std::string initialSB = configuration.getString("DMFTAB.initialSB");
        math::CFunction SA;
        math::CFunction SB;
        SA.readBinary(initialSA);
        SB.readBinary(initialSB);
        setSelfEnergyA(SA);
        setSelfEnergyB(SB);
        double scaling = -SA.getArgument(0);
        // fail if they contradict
        if(std::fabs(SA.getArgument(0) - SB.getArgument(0)) > 1e-10)
            throw std::exception();
        setScaling(scaling);
    } catch( libconfig::SettingNotFoundException& e ) {}
}

void DMFTAB::calculateGFAInt() {
    math::IntegralAB gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFTAB::gAReal, &DMFTAB::gAImag);
    gfInt.setParameterFunctions(m_selfEnergyA, m_selfEnergyB);
    gfInt.setIntegrationParameters(this);
    gfInt(m_greenFunctionA);
}

void DMFTAB::calculateGFBInt() {
    math::IntegralAB gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFTAB::gBReal, &DMFTAB::gBImag);
    gfInt.setParameterFunctions(m_selfEnergyA, m_selfEnergyB);
    gfInt.setIntegrationParameters(this);
    gfInt(m_greenFunctionB);
}

void DMFTAB::calculateFFAInt() {
    math::IntegralAB gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFTAB::fAReal, &DMFTAB::fAImag);
    gfInt.setParameterFunctions(m_selfEnergyA, m_selfEnergyB);
    gfInt.setIntegrationParameters(this);
    gfInt(m_fFunctionA);
}

void DMFTAB::calculateFFBInt() {
    math::IntegralAB gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFTAB::fBReal, &DMFTAB::fBImag);
    gfInt.setParameterFunctions(m_selfEnergyA, m_selfEnergyB);
    gfInt.setIntegrationParameters(this);
    gfInt(m_fFunctionB);
}

void DMFTAB::solve(std::string outputDirectory) {
    if(m_startWithDelta) {
        if((m_effectiveMediumA.getSize()==0) || (m_effectiveMediumB.getSize()==0)) {
            std::cerr << "no initial hybridization functions set." << std::endl;
            throw std::exception();
        }

    } else {
        if((m_selfEnergyA.getSize()==0) || (m_selfEnergyB.getSize()==0)) {
            std::cerr << "no initial self energies set." << std::endl;
            throw std::exception();
        }
    }

    double similarity = 0;
    int iteration = 1;
    while((similarity < m_tolerance) && (iteration < m_maxIterations))
    {
        std::cout << std::endl << "// ******************** \\\\" << std::endl;
        std::cout << "|| *** Iteration: " << iteration << " *** ||" << std::endl;
        std::cout << "\\\\ ******************** //" << std::endl << std::endl;
        if((iteration != 0) || (!m_startWithDelta)) {
            // use the inital self energy to calculate a lattice greens function
            calculateGFAInt();
            calculateGFBInt();
            //calculateFFAInt();
            //calculateFFBInt();
            m_greenFunctionA.write(outputDirectory + "gA_"+ utils::toString(iteration) +".dat");
            m_greenFunctionB.write(outputDirectory + "gB_"+ utils::toString(iteration) +".dat");
            //m_fFunctionA.write(outputDirectory + "fA_"+ utils::toString(iteration) +".dat");
            //m_fFunctionB.write(outputDirectory + "fB_"+ utils::toString(iteration) +".dat");

            m_effectiveMediumA.resize(m_greenFunctionA.getSize());
            m_effectiveMediumB.resize(m_greenFunctionB.getSize());

            // extract the effective medium -- trick -- A
            for(int n = 0; n < m_greenFunctionA.getSize(); n++) {
                double omega = m_greenFunctionA.getArgument(n);
                m_effectiveMediumA.set(n, omega, -(omega + I * m_delta) - m_U/2.0 - m_mu
                                     // + m_fFunctionA.getValue(n)/m_greenFunctionA.getValue(n)  );
                                       + m_selfEnergyA.getValue(n) + 1.0/m_greenFunctionA.getValue(n));
            }
            // extract the effective medium -- trick -- B
            for(int n = 0; n < m_greenFunctionB.getSize(); n++) {
                double omega = m_greenFunctionB.getArgument(n);
                m_effectiveMediumB.set(n, omega, -(omega + I * m_delta) - m_U/2.0 - m_mu
                                      //+ m_fFunctionB.getValue(n)/m_greenFunctionB.getValue(n)  );
                                       + m_selfEnergyB.getValue(n) + 1.0/m_greenFunctionB.getValue(n));
            }


            m_effectiveMediumA.rescale(1.0 / m_scaling);
            m_effectiveMediumA.write(outputDirectory + "dA_" + utils::toString(iteration) + ".dat");
            m_effectiveMediumB.rescale(1.0 / m_scaling);
            m_effectiveMediumB.write(outputDirectory + "dB_" + utils::toString(iteration) + ".dat");
        } else
            std::cout << "Skipping effective medium calculation" << std::endl;

        // feed that to the nrg as a hybridization function
        // solve impurity problem
        // update self energy

        nrg::broadening::FFTBroadener broadener;
        if(m_configuration)
            broadener.configure(*m_configuration);

        nrg::chain::HybridizationProvider chainProvider;
        if(m_configuration)
            chainProvider.configure(*m_configuration);
        chainProvider.setHybridization(m_effectiveMediumA, m_effectiveMediumB);

        //careful with that axe, Eugene.
//        if(std::abs(m_mu) < 1e-20)
//            chainProvider.setPHsymmetric(true);
//        else
//            chainProvider.setPHsymmetric(false);
        chainProvider.setPHsymmetric(false);
        chainProvider.setSZsymmetric(false);

        nrg::NRG nrg(chainProvider, broadener);
        if(m_configuration)
            nrg.configure(*m_configuration);
        nrg.setU(m_U/m_scaling);
        nrg.setEpsF((-m_U/(2.0) - m_mu) / m_scaling );
        nrg.setTemperature(m_T);
        nrg.init();
        nrg.showInfo();
        nrg.solve();
        m_n = nrg.getOccupation();
        m_m = nrg.getMagnetization();

        math::CFunction selfEnergyA;
        math::CFunction selfEnergyB;
        nrg.getSelfEnergy(selfEnergyA, selfEnergyB);
        selfEnergyA.rescale(m_scaling);
        selfEnergyB.rescale(m_scaling);
        if(iteration > 0)
            similarity = (math::similarity(selfEnergyA, m_selfEnergyA) + math::similarity(selfEnergyB, m_selfEnergyB))/2.0;
        std::cout << "iteration: " << iteration << " similarity: " << similarity << std::endl;
        m_selfEnergyA = selfEnergyA;
        m_selfEnergyB = selfEnergyB;

        m_selfEnergyA.write(outputDirectory + "sA_" + utils::toString(iteration) + ".dat");
        m_selfEnergyB.write(outputDirectory + "sB_" + utils::toString(iteration) + ".dat");

        iteration++;
    }
}


}
