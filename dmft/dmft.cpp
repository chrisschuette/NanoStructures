#include "dmft.h"
#include "cfcoefficients.h"
#include "../math/integral.h"

#include "../config/configuration.h"
#include "../nrg/chain/hybridizationprovider.h"
#include "../nrg/broadening/fftbroadener.h"
#include "../nrg/nrg.h"

#include "../mpi/openmpi.h"

#include "../utils/utils.h"

#include <fstream>

#define sqr(x) ((x)*(x))

namespace dmft {
DMFT::DMFT()
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

void DMFT::configure(config::Configuration& configuration) {
    m_configuration = &configuration;
    try {
        m_U = configuration.getDouble("DMFT.U");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_mu = configuration.getDouble("DMFT.mu");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_T = configuration.getDouble("DMFT.temperature");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_delta = configuration.getDouble("DMFT.delta");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_scaling = configuration.getDouble("DMFT.scaling");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_tolerance = configuration.getDouble("DMFT.tolerance");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_maxIterations = configuration.getInteger("DMFT.maxIterations");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_startWithDelta = configuration.getBool("DMFT.startWithDelta");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        std::string initialD = configuration.getString("DMFT.initialD");
        math::CFunction D;
        D.readBinary(initialD);
        setHybridizationFunction(D);
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        std::string initialS = configuration.getString("DMFT.initialS");
        math::CFunction S;
        S.readBinary(initialS);
        setSelfEnergy(S);
        double scaling = -S.getArgument(0);
        setScaling(scaling);
    } catch( libconfig::SettingNotFoundException& e ) {}
}


void DMFT::calculateGFInt() {
    math::Integral gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFT::gReal, &DMFT::gImag);
    gfInt.setParameterFunction(m_selfEnergy);
    gfInt.setIntegrationParameters(this);
    gfInt(m_greenFunction);
}

void DMFT::calculateFFInt() {
    math::Integral gfInt;
    gfInt.setAccuracy(1e-7,1e-7);
    gfInt.setBounds(4,-6.0,-2.0,2.0,6.0);
    gfInt.setIntegralKernels(&DMFT::fReal, &DMFT::fImag);
    gfInt.setParameterFunction(m_selfEnergy);
    gfInt.setIntegrationParameters(this);
    gfInt(m_fFunction);
}


void DMFT::calculateGFCF() {
    m_greenFunction.resize(m_selfEnergy.getSize());
    for( int i=0; i < m_selfEnergy.getSize(); i++)
    {
        double omega = m_selfEnergy.getArgument(i);
        double SelfR = m_selfEnergy.getValue(i).real();
        double SelfI = m_selfEnergy.getValue(i).imag();

        std::complex<double> Z(omega-SelfR+m_U/2.0+m_mu,-SelfI+m_delta);

        int kMAX = 301;

        double bsqr_m = bsqr[kMAX-2];

        /*         zunaechst das Abschlussglied    */
        double Re = Re_sqr(Z.real(),Z.imag());
        double Im = Im_sqr(Z.real(),Z.imag());
        Re -= 4.0 * (bsqr_m);
        double Re_neu =  Re_sqrt(Re,Im);
        double Im_neu =  Im_sqrt(Re,Im);
        Re = (Re_neu + Z.real())/2.0;
        Im = (Im_neu + Z.imag())/2.0;

        /*          Einbau in den Bruch   */

        Re_neu = bsqr[kMAX-2] * Re_calc((Z.real() - Re),Z.imag() - Im);
        Im_neu = bsqr[kMAX-2] * Im_calc((Z.real() - Re),Z.imag() - Im);
        Re = Re_neu;
        Im = Im_neu;

        for (int j=(kMAX-3);j>=0;j--) {
            Re_neu = bsqr[j] * Re_calc((Z.real()  - Re),(Z.imag() - Im));
            Im_neu = bsqr[j] * Im_calc((Z.real()  - Re),(Z.imag() - Im));

            Re = Re_neu;
            Im = Im_neu;
        }

        Re_neu =  Re_calc((Z.real() - Re),(Z.imag() - Im));
        Im_neu =  Im_calc((Z.real() - Re),(Z.imag() - Im));

        double norm = 1.0;

        Re = norm * Re_neu / M_PI;
        Im = norm * Im_neu / M_PI;

        std::complex<double> G(M_PI*Re, M_PI*Im);

        m_greenFunction.set(i, omega, G);
    }
}

void DMFT::solve(std::string outputDirectory) {
/*    if((m_selfEnergy.getSize()==0) && (m_effectiveMedium.getSize() == 0))
    {
        std::cerr << "no initial self energy nor hybridization function set." << std::endl;
        throw std::exception();
    } else {
        std::cout << m_selfEnergy.getSize() << std::endl;
        std::cout << m_effectiveMedium.getSize() << std::endl;
    }*/
    double n;
    if(m_startWithDelta) {
        if(m_effectiveMedium.getSize()==0) {
            std::cerr << "no initial hybridization function set." << std::endl;
            throw std::exception();
        }

    } else {
        if(m_selfEnergy.getSize()==0) {
            std::cerr << "no initial self energy set." << std::endl;
            throw std::exception();
        }
    }


    m_selfEnergy.write("sInitial.dat");
    double similarity = 0;
    int iteration = 1;
    while((similarity < m_tolerance) && (iteration < m_maxIterations))
    {

        std::cout << std::endl << "// ******************** \\\\" << std::endl;
        std::cout << "|| *** Iteration: " << iteration << " *** ||" << std::endl;
        std::cout << "\\\\ ******************** //" << std::endl << std::endl;

        if((iteration != 0) || (!m_startWithDelta)) {
            // use the inital self energy to calculate a lattice greens function
            calculateGFInt();
            calculateFFInt();
            m_greenFunction.write(outputDirectory + "g_"+ utils::toString(iteration) +".dat");
            m_fFunction.write(outputDirectory + "f_"+ utils::toString(iteration) +".dat");

            m_effectiveMedium.resize(m_greenFunction.getSize());

            // extract the effective medium -- trick
            for(int n = 0; n < m_greenFunction.getSize(); n++) {
                double omega = m_greenFunction.getArgument(n);
                m_effectiveMedium.set(n, omega, -(omega + I * m_delta) - m_U/2.0 - m_mu
                                      + m_fFunction.getValue(n)/m_greenFunction.getValue(n)  );
            }
            m_effectiveMedium.rescale(1.0 / m_scaling);
            m_effectiveMedium.write(outputDirectory + "d_" + utils::toString(iteration) + ".dat");
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
        chainProvider.setHybridization(m_effectiveMedium, m_effectiveMedium);

        //careful with that axe, Eugene.
        if(std::abs(m_mu) < 1e-20)
            chainProvider.setPHsymmetric(true);
        else
            chainProvider.setPHsymmetric(false);
        chainProvider.setSZsymmetric(true);

        nrg::NRG nrg(chainProvider, broadener);
        if(m_configuration)
            nrg.configure(*m_configuration);
        nrg.setU(m_U/m_scaling);
        nrg.setEpsF((-m_U/(2.0) - m_mu) / m_scaling );
        nrg.setTemperature(m_T);
        nrg.init();
        nrg.showInfo();
        nrg.solve();
        n = nrg.getOccupation();

        math::CFunction selfEnergy;
        nrg.getSelfEnergy(selfEnergy);
        selfEnergy.rescale(m_scaling);
        if(iteration > 0)
            similarity = math::similarity(selfEnergy, m_selfEnergy);
        std::cout << "iteration: " << iteration << " similarity: " << similarity << std::endl;
        m_selfEnergy = selfEnergy;

        m_selfEnergy.write(outputDirectory + "s_" + utils::toString(iteration) + ".dat");
        iteration++;
    }
    if(iteration > m_maxIterations) {
        std::cerr << "**** DID NOT CONVERGE ****" << std::endl;
    }
    std::ofstream nFile(std::string("n.dat").c_str(), std::ios::app);
    nFile << m_mu << " " << n << std::endl;
    nFile.close();
}


double DMFT::Re_calc(double Re, double Im)
{
    double erg;

    erg = Re/(sqr(Re) + sqr(Im));
    return erg;
}

double DMFT::Im_calc(double Re, double Im)
{
    double erg;

    erg = -1.0 * Im/(sqr(Re) + sqr(Im));
    return erg;
}

double DMFT::Re_sqrt(double Re, double Im)
{
    double erg,betrag,phase,Pi;
    Pi = acos(-1.0);

    betrag = sqrt(sqr(Re) + sqr(Im));
    phase = atan2(Im,Re);
    if (phase < 0) {
        phase = 2.0 * Pi + phase;
    }
    erg =  -1.0 * sqrt(betrag) * cos(phase/2.0);
    return erg;
}


double DMFT::Im_sqrt(double Re, double Im)
{
    double erg,betrag,phase,Pi;
    Pi = acos(-1.0);

    betrag = sqrt(sqr(Re) + sqr(Im));
    phase = atan2(Im,Re);
    if (phase < 0) {
        phase = 2.0 * Pi + phase;
    }
    erg = -1.0 * sqrt(betrag) * sin(phase/2.0);
    return erg;
}

double DMFT::Re_sqr(double Re, double Im)
{
    double erg;

    erg = sqr(Re) - sqr(Im);
    return erg;
}

double DMFT::Im_sqr(double Re, double Im)
{
    double erg;

    erg = 2.0 * Re * Im;
    return erg;
}

DMFT::~DMFT()
{

}

}
