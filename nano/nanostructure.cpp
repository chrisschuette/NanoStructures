#include "nanostructure.h"
#include "../math/integral.h"
#include "../utils/utils.h"
#include "../utils/utils_io.h"
#include "../config/configuration.h"
#include "../nrg/broadening/fftbroadener.h"
#include "../nrg/chain/hybridizationprovider.h"
#include "../nrg/nrg.h"
#include "../mpi/openmpi.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <libconfig.h++>

namespace nano {
NanoStructure::NanoStructure()
    : m_nLayers(0)
    , m_n(0)
    , m_selfEnergies(0)
    , m_greensFunctions(0)
    , m_fFunctions(0)
    , m_effectiveMedia(0)
    , m_U(0)
    , m_mu(0)
    , m_T(1e-8)
    , m_doECR(false)
    , m_rho(0)
    , m_rhoBackground(0)
    , m_V(0)
    , m_symmetricNano(false)
    , m_tL(0)
    , m_tR(0)
    , m_similarities()
    , m_delta(1e-2)
    , m_scaling(100.0)
    , m_alphaV(0.05)
    , m_tolerance(0.995)
    , m_maxIterations(20)
    , m_L(0)
    , m_dL(0)
    , m_R(0)
    , m_dR(0)
    , m_temp(0)
    , m_configuration(0)
    , m_calculateRhoXX(false)
    , m_calculateRhoXY(false)
    , m_verbose(false)
    , m_iteration(1)
{
}

NanoStructure::NanoStructure(int n)
    : m_nLayers(n)
    , m_T(1e-8)
    , m_doECR(false)
    , m_symmetricNano(false)
    , m_delta(1e-2)
    , m_scaling(100.0)
    , m_alphaV(0.05)
    , m_tolerance(0.995)
    , m_maxIterations(20)
    , m_configuration(0)
    , m_calculateRhoXX(false)
    , m_calculateRhoXY(false)
    , m_verbose(false)
    , m_iteration(1)
{
    allocate(n);
    // initialize with default values
    for(int i = 0; i < n+2; i++) {
        m_U[i] = 0.0;
        m_mu[i] = 0.0;
        m_rho[i] = 1.0;
        m_rhoBackground[i] = 1.0;
        m_V[i] = 0.0;
        m_tL[i] = 1.0;
        m_tR[i] = 1.0;
        m_similarities[i] = 0.0;
    }
}

NanoStructure::NanoStructure(std::string filename)
    : m_nLayers(0)
    , m_selfEnergies(0)
    , m_greensFunctions(0)
    , m_fFunctions(0)
    , m_effectiveMedia(0)
    , m_U(0)
    , m_mu(0)
    , m_T(1e-8)
    , m_doECR(false)
    , m_rho(0)
    , m_rhoBackground(0)
    , m_V(0)
    , m_symmetricNano(false)
    , m_tL(0)
    , m_tR(0)
    , m_similarities(0)
    , m_delta(1e-2)
    , m_scaling(100.0)
    , m_alphaV(0.05)
    , m_tolerance(0.995)
    , m_maxIterations(20)
    , m_L(0)
    , m_dL(0)
    , m_R(0)
    , m_dR(0)
    , m_temp(0)
    , m_configuration(0)
    , m_calculateRhoXX(false)
    , m_calculateRhoXY(false)
    , m_verbose(false)
    , m_iteration(1)
{
    // step 1 -- read the structure file
    libconfig::Config structureFile;
    structureFile.readFile(filename.c_str());
    tPMap sections = parseStructure(structureFile);

    // look up structure section
    libconfig::Setting& structure = structureFile.lookup("structure");
    libconfig::Setting& layers = structureFile.lookup("structure.layers");
    assert(structure.isGroup());
    assert(layers.isGroup());

    int N = 0; // two leads

    // iterate over sections and calculate the total number of layers
    for(int i = 0; i < layers.getLength(); i++) {
        assert(layers[i].isList() && layers[i].getLength() == 2);
        N += (int) layers[i][1];
    }
    std::cout << N << " layers." << std::endl;

    // allocate memory
    m_nLayers = N;
    allocate(N);

    // left lead
    std::string L = structureFile.lookup("structure.L");
    tPMap::const_iterator i = sections.find(L);
    if(i == sections.end()) {
        std::cerr << "Could not find left lead." << std::endl;
        throw std::exception();
    }
    std::cout << "L = " << L << std::endl;
    m_V[0] = (*i).second.initialV;
    m_U[0] = (*i).second.U;
    m_mu[0] = (*i).second.mu;
    m_rhoBackground[0] = (*i).second.rhoB;
    m_rho[0] = (*i).second.rhoB;
    if((*i).second.selfEnergyFile != "")
        m_selfEnergies[0].readBinary((*i).second.selfEnergyFile);
    if((*i).second.occupancyFile != "")
        m_n[0].read((*i).second.occupancyFile);

    //iterate over sections
    int currentLayer = 1;
    for(int j = 0; j < layers.getLength(); j++) {
        std::string material = layers[j][0];
        int quantity = layers[j][1];
        assert(quantity >= 0);
        i = sections.find(material);
        if(i == sections.end()) {
            std::cerr << "Could not find section " << material << "." << std::endl;
            throw std::exception();
        }
        std::cout << quantity << "x " << material << std::endl;
        for(int c = 0; c < quantity; c++) {
            m_V[currentLayer] = (*i).second.initialV;
            m_U[currentLayer] = (*i).second.U;
            m_mu[currentLayer] = (*i).second.mu;
            m_rhoBackground[currentLayer] = (*i).second.rhoB;
            m_rho[currentLayer] = (*i).second.rhoB;
            if((*i).second.selfEnergyFile != "")
                m_selfEnergies[currentLayer].readBinary((*i).second.selfEnergyFile);
            if((*i).second.occupancyFile != "")
                m_n[currentLayer].read((*i).second.occupancyFile);
            currentLayer++;
        }
    }

    // right lead
    std::string R = structureFile.lookup("structure.R");
    i = sections.find(R);
    if(i == sections.end()) {
        std::cerr << "Could not find right lead." << std::endl;
        throw std::exception();
    }
    std::cout << "R = " << R << std::endl;
    m_V[N+1] = (*i).second.initialV;
    m_U[N+1] = (*i).second.U;
    m_mu[N+1] = (*i).second.mu;
    m_rhoBackground[N+1] = (*i).second.rhoB;
    m_rho[N+1] = (*i).second.rhoB;
    if((*i).second.selfEnergyFile != "")
        m_selfEnergies[N+1].readBinary((*i).second.selfEnergyFile);
    if((*i).second.occupancyFile != "")
        m_n[N+1].read((*i).second.occupancyFile);

    // initialize with default values
    for(int i = 0; i < N+2; i++) {
        m_tL[i] = 1.0;
        m_tR[i] = 1.0;
        m_similarities[i] = 0.0;
    }
}

NanoStructure::tPMap NanoStructure::parseStructure(const libconfig::Config& structure) {
    std::map<std::string, nano::LayerParameters> sections;
    libconfig::Setting& root = structure.getRoot();
    for(int i = 0; i < root.getLength(); i++) {
        if(root[i].isGroup() && (std::string(root[i].getName()) != "structure")) {
            libconfig::Setting& section = root[i];
            std::cout << "Parsing section " << section.getName() << std::endl;
            nano::LayerParameters parameters;
            section.lookupValue("U", parameters.U);
            section.lookupValue("mu", parameters.mu);
            section.lookupValue("rhoB", parameters.rhoB);
            section.lookupValue("S", parameters.selfEnergyFile);
            section.lookupValue("V", parameters.initialV);
            section.lookupValue("N", parameters.occupancyFile);
            // ...
            sections.insert(std::pair<std::string, nano::LayerParameters>(section.getName(), parameters));
        }
    }
    return sections;
}

void NanoStructure::allocate(int N) {
    m_n=new math::Function[N+2];
    m_selfEnergies=new math::CFunction[N+2];
    m_greensFunctions=new math::CFunction[N+2];
    m_fFunctions=new math::CFunction[N+2];
    m_effectiveMedia=new math::CFunction[N+2];
    m_U=new double[N+2];
    m_mu=new double[N+2];
    m_rho=new double[N+2];
    m_rhoBackground=new double[N+2];
    m_V=new double[N+2];
    m_tL=new double[N+2];
    m_tR=new double[N+2];
    m_similarities=new double[N+2];
    m_L=new std::complex<double>[N+2];
    m_dL=new std::complex<double>[N+2];
    m_R=new std::complex<double>[N+2];
    m_dR=new std::complex<double>[N+2];
    m_temp=new std::complex<double>[N+2];
}

void NanoStructure::free() {
    if(m_n)
        delete [] m_n;
    if(m_selfEnergies)
        delete [] m_selfEnergies;
    m_selfEnergies = 0;
    if(m_greensFunctions)
        delete [] m_greensFunctions;
    m_greensFunctions = 0;
    if(m_effectiveMedia)
        delete [] m_effectiveMedia;
    m_effectiveMedia = 0;
    if(m_U)
        delete [] m_U;
    m_U = 0;
    if(m_mu)
        delete [] m_mu;
    m_mu = 0;
    if(m_rho)
        delete [] m_rho;
    m_rho = 0;
    if(m_rhoBackground)
        delete [] m_rhoBackground;
    m_rhoBackground = 0;
    if(m_V)
        delete [] m_V;
    m_V = 0;
    if(m_tL)
        delete [] m_tL;
    m_tL = 0;
    if(m_tR)
        delete [] m_tR;
    m_tR = 0;
    if(m_similarities)
        delete [] m_similarities;
    m_similarities = 0;
    if(m_temp)
        delete [] m_temp;
    m_temp = 0;
    if(m_L)
        delete [] m_L;
    m_L = 0;
    if(m_dL)
        delete [] m_dL;
    m_dL = 0;
    if(m_R)
        delete [] m_R;
    m_R = 0;
    if(m_dR)
        delete [] m_dR;
    m_dR = 0;
}

NanoStructure::~NanoStructure() {
    free();
}

void NanoStructure::configure(config::Configuration& config) {
    m_configuration = &config;
    try {
      m_T = m_configuration->getDouble("NANO.T");
    } catch( libconfig::SettingNotFoundException& e ) {}
}


void NanoStructure::precalc_dL(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);
    std::complex<double> a = (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0]))/2.0;

    std::complex<double> b = m_L[0] - a;
    m_dL[0] = -0.5 - 2.0 * a / (4.0 * b);

    assert(nLayer <= m_nLayers);
    for(int n = 1; n <= nLayer; n++) {
        assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
        m_dL[n] = -1.0 + m_dL[n-1] / std::pow(m_L[n-1],2);
    }
}

void NanoStructure::precalc_dR(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);
    std::complex<double> a = (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1]))/2.0;

    std::complex<double> b = m_R[m_nLayers+1] - a;
    m_dR[m_nLayers+1] = -0.5 - 2.0 * a / (4.0 * b);

    assert(nLayer <= m_nLayers);
    for(int n = m_nLayers; n >= nLayer; n--) {
        assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
        m_dR[n] = -1.0 + m_dR[n+1] / std::pow(m_R[n+1],2);
    }
}

void NanoStructure::precalc_L(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);

    std::complex<double> a = (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0]))/2.0;
    std::complex<double> b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0])) , 2 ) - 4.0 );

    m_L[0] = a + b;
    if(m_L[0].imag() < 0)
        m_L[0] = a - b;

    assert(nLayer <= m_nLayers);
    for(int n = 1; n <= nLayer; n++) {
        assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
        std::complex<double> a = z - m_selfEnergies[n].getValue(nOmega) - (epsk - m_U[n]/2.0 - m_mu[n] + m_V[n]);
        m_L[n] = a - 1.0 / m_L[n-1];
    }
}

void NanoStructure::precalc_R(int nLayer, int nOmega, double epsk) {
    assert(m_selfEnergies[0].getSize() == m_selfEnergies[nLayer].getSize());
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);

    std::complex<double> a = (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1]))/2.0;
    std::complex<double> b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1])) , 2 ) - 4.0 );

    m_R[m_nLayers+1] = a + b;
    if(m_R[m_nLayers+1].imag() < 0)
        m_R[m_nLayers+1] = a - b;

    assert(nLayer > 0);
    for(int n = m_nLayers; n >= nLayer; n--) {
        assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
        std::complex<double> a = z - m_selfEnergies[n].getValue(nOmega) - (epsk - m_U[n]/2.0 - m_mu[n] + m_V[n]);
        m_R[n] = a - 1.0 / m_R[n+1];
    }
}

std::complex<double> NanoStructure::left(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);
    std::complex<double> L;

    std::complex<double> a = (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0]))/2.0;
    std::complex<double> b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0])) , 2 ) - 4.0 );

    L = a + b;
    if(L.imag() < 0)
        L = a - b;

    if(nLayer == 0)
        return L;
    else {
        assert(nLayer <= m_nLayers);
        for(int n = 1; n <= nLayer; n++) {
            assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
            std::complex<double> a = z - m_selfEnergies[n].getValue(nOmega) - (epsk - m_U[n]/2.0 - m_mu[n] + m_V[n]);
            L = a - 1.0 /  L;
        }
    }
    return L;
}

std::complex<double> NanoStructure::right(int nLayer, int nOmega, double epsk) {
    assert(m_selfEnergies[0].getSize() == m_selfEnergies[nLayer].getSize());
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);

    std::complex<double> R;

    std::complex<double> a = (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1]))/2.0;
    std::complex<double> b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1])) , 2 ) - 4.0 );

    R = a + b;
    if(R.imag() < 0)
        R = a - b;

    if(nLayer == (m_nLayers + 1))
        return R;
    else {
        assert(nLayer > 0);
        for(int n = m_nLayers; n >= nLayer; n--) {
            assert(m_selfEnergies[n].getSize() == m_selfEnergies[nLayer].getSize());
            std::complex<double> a = z - m_selfEnergies[n].getValue(nOmega) - (epsk - m_U[n]/2.0 - m_mu[n] + m_V[n]);
            R = a - 1.0 / R;
        }
    }
    return R;
}

double NanoStructure::Pi(int alpha, int gamma, int delta, int nOmega, double epsk) {
    int max  = std::max(gamma, delta);
    int min = std::min(gamma, delta);
    int maxIndex = std::max(alpha , max);
    int minIndex = std::min(alpha , min);

    //precalculate stuff
    precalc_L(maxIndex, nOmega, epsk);
    precalc_dL(max, nOmega, epsk);
    precalc_R(minIndex, nOmega, epsk);
    precalc_dR(min, nOmega, epsk);

    // figure out G'(gamma, delta)

    // 1) calculate G'(max(gamma,delta), max(gamma,delta))
    std::complex<double> z(m_selfEnergies[max].getArgument(nOmega), m_delta);
    std::complex<double> a = z - m_selfEnergies[max].getValue(nOmega) - (epsk - m_U[max]/2.0 - m_mu[max] + m_V[max]);
    std::complex<double> G_max_max = 1.0 / (m_L[max] + m_R[max] - a);
    std::complex<double> dG_max_max = - (1.0 + m_dL[max] + m_dR[max]) * std::pow(G_max_max, 2);

    // 2) use the L functions and their derivatives to move it back!
    std::complex<double> G_min_max = 0.0;
    std::complex<double> dG_min_max = 0.0;
    if(max != min) {
        G_min_max = G_max_max;
        for(int i = max-1; i>= min; i--)
            G_min_max *= -(1.0 / m_L[i]);

        // G' term
        std::complex<double> contribution = dG_max_max;
        for(int i = max-1; i>= min; i--)
            contribution *= -(1.0 / m_L[i]);
        dG_min_max += contribution;
        // the others
        for(int i = max-1; i>= min; i--) {
            contribution = G_max_max;
            for(int j = max-1; j > i; j--) {
                contribution *= -(1.0 / m_L[i]);
            }
            contribution *= m_dL[i] / std::pow(m_L[i],2);
            for(int j = i-1; j >= min; j--) {
                contribution *= -(1.0 / m_L[i]);
            }
            dG_min_max += contribution;
        }
    } else {
        dG_min_max = dG_max_max;
        G_min_max = G_max_max;
    }

    int maxAlphaGamma = std::max(alpha, gamma);
    int maxAlphaDelta = std::max(alpha, delta);
    int minAlphaGamma = std::min(alpha, gamma);
    int minAlphaDelta = std::min(alpha, delta);

    a = z - m_selfEnergies[maxAlphaGamma].getValue(nOmega) - (epsk - m_U[maxAlphaGamma]/2.0 - m_mu[maxAlphaGamma] + m_V[maxAlphaGamma]);
    std::complex<double> G_alpha_gamma = 1.0 / (m_L[maxAlphaGamma] + m_R[maxAlphaGamma] - a);
    if(alpha != gamma)
        for(int i = maxAlphaGamma-1; i>= minAlphaGamma; i--)
            G_alpha_gamma *= -(1.0 / m_L[i]);

    a = z - m_selfEnergies[maxAlphaDelta].getValue(nOmega) - (epsk - m_U[maxAlphaDelta]/2.0 - m_mu[maxAlphaDelta] + m_V[maxAlphaDelta]);
    std::complex<double> G_alpha_delta = 1.0 / (m_L[maxAlphaDelta] + m_R[maxAlphaDelta] - a);
    if(alpha != delta)
        for(int i = maxAlphaDelta-1; i>= minAlphaDelta; i--)
            G_alpha_delta *= -(1.0 / m_L[i]);

    return -1.0 / M_PI * ( std::imag(G_alpha_gamma * G_alpha_delta * dG_min_max) * math::physics::fermi(1/m_T, epsk)
                           + std::real(G_alpha_gamma * G_alpha_delta) * std::imag(G_min_max) * math::physics::dfermi(1/m_T, epsk)
                           );
}

std::complex<double> NanoStructure::G(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);
    std::complex<double> a = (z - m_selfEnergies[nLayer].getValue(nOmega) - (epsk - m_U[nLayer]/2.0 - m_mu[nLayer] + m_V[nLayer]));
    return 1.0 / (left(nLayer, nOmega, epsk) + right(nLayer, nOmega, epsk) - a);
}

std::complex<double> NanoStructure::F(int nLayer, int nOmega, double epsk) {
    std::complex<double> z(m_selfEnergies[nLayer].getArgument(nOmega), m_delta);
    std::complex<double> a = (z - m_selfEnergies[nLayer].getValue(nOmega) - (epsk - m_U[nLayer]/2.0 - m_mu[nLayer] + m_V[nLayer]));
    std::complex<double> G = 1.0 / (left(nLayer, nOmega, epsk) + right(nLayer, nOmega, epsk) - a);
    return m_selfEnergies[nLayer].getValue(nOmega) * G + 1.0;
}

std::complex<double> NanoStructure::GOff(int alpha, int beta, int nOmega, double epsk) {
    if(alpha == beta)
        return G(alpha, nOmega, epsk);
    else if(alpha > beta)
        utils::swap(alpha, beta);

    std::complex<double> z(m_selfEnergies[beta].getArgument(nOmega), m_delta);

    // calculate L_-inf and R_inf
    /*
    //calculate right slab
    std::complex<double> a;
    a = (z - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0)) / 2.0;
    std::complex<double> b;
    b = std::sqrt(std::pow((z - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 )), 2.0) - 4 * m_tL[m_nLayers+1] * m_tR[m_nLayers+1]) / 2.0;
    double sign;
    std::complex<double> o1 = a + b;
    if (o1.imag() < 0)
        sign = -1.0;
    else
        sign = 1.0;
    temp[m_nLayers+1] = a + sign * b;
    */
    std::complex<double> a = (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1]))/2.0;
    std::complex<double> b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[m_nLayers+1].getValue(nOmega) - (epsk - m_U[m_nLayers+1]/2.0 - m_mu[m_nLayers+1] + m_V[m_nLayers+1])) , 2 ) - 4.0 );

    m_temp[m_nLayers+1] = a + b;
    if(m_temp[m_nLayers+1].imag() < 0)
        m_temp[m_nLayers+1] = a - b;
    /*
    //calculate left slab
    a = (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0)) / 2.0;
    b = std::sqrt(std::pow((z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0)), 2.0) - 4 * m_tR[0] * m_tR[0]) / 2.0;
    o1 = a + b;
    if (o1.imag() < 0)
        sign = -1.0;
    else
        sign = 1.0;
    temp[0] = a + sign * b;
    */
    a = (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0 - m_mu[0] + m_V[0]))/2.0;
    b = 0.5 * std::sqrt( std::pow( (z  - m_selfEnergies[0].getValue(nOmega) - (epsk - m_U[0]/2.0) - m_mu[0] + m_V[0]) , 2 ) - 4.0 );

    m_temp[0] = a + b;
    if(m_temp[0].imag() < 0)
        m_temp[0] = a - b;

    // calculate all left functions < beta
    for(int i = 1; i < beta; i++) {
        a = z - m_selfEnergies[i].getValue(nOmega) - (epsk - m_U[i]/2.0 - m_mu[i] + m_V[i]);
        if (m_tL[i] * m_tR[i-1] == 0.0)
            b = std::complex<double>(0.0, 0.0);
        else
            b = m_tL[i] * m_tR[i-1] / m_temp[i-1];
        m_temp[i] = a - b;
    }

    //calculate all right functions > beta
    for(int i = m_nLayers; i > beta; i--) {
        a = z - m_selfEnergies[i].getValue(nOmega) - (epsk - m_U[i]/2.0 - m_mu[i] + m_V[i]);
        if (m_tR[i] * m_tL[i+1] == 0.0)
            b = std::complex<double>(0.0, 0.0);
        else
            b = m_tR[i] * m_tL[i+1] / m_temp[i+1];
        m_temp[i] = a - b;
    }

    // calculate G_beta,beta

    //left
    a = z - m_selfEnergies[beta].getValue(nOmega) - (epsk - m_U[beta]/2.0 - m_mu[beta] + m_V[beta]);
    if (m_tL[beta] * m_tR[beta-1] == 0.0)
        b = std::complex<double>(0.0, 0.0);
    else
        b = m_tL[beta] * m_tR[beta-1] / m_temp[beta-1];
    std::complex<double> L = a - b;

    // right
    if (m_tR[beta] * m_tL[beta+1] == 0.0)
        b = std::complex<double>(0.0, 0.0);
    else
        b = m_tR[beta] * m_tL[beta+1] / m_temp[beta+1];
    std::complex<double> R = a - b;

    // G_beta,beta
    m_temp[beta] = 1.0 / (L + R - a);

    // calculate G_alpha,beta
    std::complex<double> G = m_temp[beta];
    for(int i = beta-1; i>= alpha; i--) {
        G *= -(1.0 / m_temp[i]);
    }
    return G;
}

void NanoStructure::F_local(int nLayer, math::CFunction& fLocal, mpi::OpenMPI& mpi) {
    IntegrationParametersG parameters;
    parameters.nLayer = nLayer;
    parameters.structure = this;
    math::Integral gInt;
    gInt.setAccuracy(0.0,1e-3);
    gInt.setBounds(5,-4.0,-2.0,0.0,2.0,4.0);
    gInt.setIntegralKernels(&NanoStructure::fReal, &NanoStructure::fImag);
    gInt.setIntegrationParameters(&parameters);
    gInt.setParameterFunction(m_selfEnergies[nLayer]);
    gInt(fLocal, mpi);
}

void NanoStructure::G_local(int nLayer, math::CFunction& gLocal, mpi::OpenMPI& mpi) {
    IntegrationParametersG parameters;
    parameters.nLayer = nLayer;
    parameters.structure = this;
    math::Integral gInt;
    gInt.setAccuracy(0.0,1e-3);
    gInt.setBounds(5,-4.0,-2.0,0.0,2.0,4.0);
    gInt.setIntegralKernels(&NanoStructure::gReal, &NanoStructure::gImag);
    gInt.setIntegrationParameters(&parameters);
    gInt.setParameterFunction(m_selfEnergies[nLayer]);
    gInt(gLocal, mpi);
}

void NanoStructure::F_local(int nLayer, math::CFunction& fLocal) {
    IntegrationParametersG parameters;
    parameters.nLayer = nLayer;
    parameters.structure = this;
    math::Integral gInt;
    gInt.setAccuracy(0.0,1e-3);
    gInt.setBounds(5,-4.0,-2.0,0.0,2.0,4.0);
    gInt.setIntegralKernels(&NanoStructure::fReal, &NanoStructure::fImag);
    gInt.setIntegrationParameters(&parameters);
    gInt.setParameterFunction(m_selfEnergies[nLayer]);
    gInt(fLocal);
}

void NanoStructure::G_local(int nLayer, math::CFunction& gLocal) {
    IntegrationParametersG parameters;
    parameters.nLayer = nLayer;
    parameters.structure = this;
    math::Integral gInt;
    gInt.setAccuracy(0.0,1e-3);
    gInt.setBounds(5,-4.0,-2.0,0.0,2.0,4.0);
    gInt.setIntegralKernels(&NanoStructure::gReal, &NanoStructure::gImag);
    gInt.setIntegrationParameters(&parameters);
    gInt.setParameterFunction(m_selfEnergies[nLayer]);
    gInt(gLocal);
}


void NanoStructure::setSelfEnergies(const math::CFunction& selfEnergy) {
    for(int n = 0; n < m_nLayers+2; n++)
        m_selfEnergies[n] = selfEnergy;
}

void NanoStructure::setHoppings(double tL, double tR) {
    for(int n = 0; n < m_nLayers+2; n++) {
        m_tL[n] = tL;
        m_tR[n] = tR;
    }
}

void NanoStructure::setUs(double U) {
    for(int n = 0; n < m_nLayers+2; n++) {
        m_U[n] = U;
    }
}

bool NanoStructure::isConverged() {
    for(int n = 1; n <= m_nLayers; n++)
        if(m_similarities[n] < m_tolerance)
            return false;
    if(m_doECR) {
        double tECpL = std::abs(totalExcessCharge(0.0)/((double) m_nLayers));
        if(tECpL > 1e-6) {
            std::cout << "total Excess charge avoids convergence: " << tECpL << std::endl;
            return false;
        }
    }
    return true;
}

double NanoStructure::doNRG(int n) {
    nrg::broadening::FFTBroadener broadener;
    if(m_configuration != 0)
        broadener.configure(*m_configuration);
    nrg::chain::HybridizationProvider chainProvider;
    if(m_configuration != 0)
        chainProvider.configure(*m_configuration);
    chainProvider.setHybridization(m_effectiveMedia[n], m_effectiveMedia[n]);
    chainProvider.setPHsymmetric(false);
    chainProvider.setSZsymmetric(true);

    nrg::NRG nrg(chainProvider, broadener);
    if(m_configuration != 0)
        nrg.configure(*m_configuration);
    nrg.setU(m_U[n]/m_scaling);
    nrg.setEpsF((-m_U[n]/2.0 - m_mu[n] + m_V[n])/m_scaling);
    nrg.setTemperature(m_T); // do I have to rescale the temperature too?
    nrg.init();
    if(m_verbose)
        nrg.showInfo();
    nrg.solve(true);

    math::CFunction selfEnergy;
    nrg.getSelfEnergy(selfEnergy);
    selfEnergy.rescale(m_scaling);
    m_selfEnergies[n] = selfEnergy;
    return nrg.getOccupation();
}

void NanoStructure::calculateEffectiveMedia(mpi::OpenMPI& mpi) {
    // calculate local Greensfunctions and occupations
    for(int n = 1; n <= m_nLayers; n++) {
        HEAD std::cout << "Calculating Green's function for layer " << n << std::endl;
        math::CFunction F;
        math::CFunction G;

        // calculat F(w) and G(w)
        F_local(n, F, mpi);
        G_local(n, G, mpi);

        //update similarity record
        if(m_iteration >= 2)
            m_similarities[n] = math::similarity(G, m_greensFunctions[n]);
        HEAD std::cout << "Layer " << n << " similarity: " << m_similarities[n] << " (" << m_tolerance << ")" << std::endl;
        m_greensFunctions[n] = G;
        m_fFunctions[n] = F;
        //HEAD m_greensFunctions[n].write("g_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");
    }

    //calculate effective media
    for(int n = 1; n <= m_nLayers; n++) {
        HEAD std::cout << "Calculating effective medium for layer " << n << std::endl;

        m_effectiveMedia[n].resize(m_greensFunctions[n].getSize());
        for(int nOmega = 0; nOmega < m_greensFunctions[m_nLayers].getSize(); nOmega++) {
            double omega = m_greensFunctions[n].getArgument(nOmega);
            m_effectiveMedia[n].set(nOmega, omega, -(omega + std::complex<double>(0.0, 1.0) * m_delta) - m_U[n]/2.0 - m_mu[n] + m_V[n]
                                    + m_fFunctions[n].getValue(nOmega) / m_greensFunctions[n].getValue(nOmega) );
        }
        m_effectiveMedia[n].rescale(1.0 / m_scaling);
        //HEAD m_effectiveMedia[n].write("d_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");
    }
}

void NanoStructure::solve(mpi::OpenMPI& mpi) {
    HEAD std::cout << " ** Solving NanoStructure using " << mpi.getSize() << " node" << (mpi.getSize()==1 ? "" : "s" ) << ". ** " << std::endl;
    HEAD std::cout << "ECR: " << (m_doECR ? "true" : "false") << std::endl;

    // solve the poisson problem
    if(m_doECR)
        initializePotentials();

    while( (m_iteration <= m_maxIterations) ) {

        HEAD std::cout << "Iteration: " << m_iteration << std::endl;

        HEAD {
            if(m_iteration % 50 == 0)
                for(int i = 1; i <= m_nLayers; i++) {
                    m_greensFunctions[i].write("g_" + utils::toString(i) + "_" + utils::toString(m_iteration) + ".txt");
                    m_selfEnergies[i].write("s_" + utils::toString(i) + "_" + utils::toString(m_iteration) + ".txt");
                }
        }

        // step 0 -- perform consistency checks
        if(!areAllSelfEnergySizesEqual() || !isScalingCompatibleWithSelfEnergies()) {
            std::cerr << "Consistency violated." << std::endl;
            throw std::exception();
        }

        // step 1 -- update the Hartree potentials (optional)
        // * no need to parallelize this.
        if(m_doECR)
            updatePotentials(mpi);

        // step 2 -- calculate the effective media for the impurity problems
        // * parallelize the integration
        calculateEffectiveMedia(mpi);

        // step 2b -- if the heterostruture is symmetric under PH + z->Nz-z average
        if(m_symmetricNano)
            for(int n = 1; n <= (m_nLayers/2 + m_nLayers%2); n++) {
                math::CFunction avg = (m_effectiveMedia[n] + m_effectiveMedia[m_nLayers-n+1].transformPH())/2.0;
                m_effectiveMedia[n] = avg;
                m_effectiveMedia[m_nLayers-n+1] = avg.transformPH();
            }

        //step 3 -- solve the impurity problem
        // * we will divide the problems among the processors in a round robin fashion
        for(int n = mpi.getRank()+1; n <= m_nLayers; n += mpi.getSize())
            m_rho[n] = doNRG(n);

        mpi.sync();

        // step 4 -- exchange information ( self energies and occupations )
        for(int n = 1; n <= m_nLayers; n++) {
            int source = (n - 1) % mpi.getSize();
            mpi.sync(m_selfEnergies[n], source);
            mpi.sync(m_rho[n], source);
        }

        HEAD writeBinary("structure");

        m_iteration++;
    }

    if(m_doECR)
        HEAD std::cout << "Final excess charge: " << totalExcessCharge(0.0) << std::endl;

    if(m_iteration > m_maxIterations)
        HEAD std::cout << "WARNING: Calculation did not converge." << std::endl;

    // output
    HEAD for(int n = 0; n < m_nLayers+2; n++) {
        m_greensFunctions[n].write("G_" + utils::toString(n) + ".txt");
        m_selfEnergies[n].write("S_" + utils::toString(n) + ".txt");
    }

    //post processing
    if(m_calculateRhoXX) {
        double sigmaXX[m_nLayers+2];
        double piXX[m_nLayers+2][m_nLayers+2];

        int counter = 1;
        for(int alpha = 1; alpha <= m_nLayers; alpha++)
            for(int beta = alpha; beta <= m_nLayers; beta++) {
                piXX[alpha][beta] = calculateSigmaXX(alpha, beta, mpi);
                HEAD std::cout << alpha << " " << beta << " " << piXX[alpha][beta] << " " << counter << "/" << m_nLayers * (m_nLayers+1) / 2 << std::endl;
                counter++;
            }

	std::ofstream sigmaFile;
	std::ofstream sigmaMatrix;
	HEAD sigmaFile.open("sigmaXX.dat");
	HEAD sigmaMatrix.open("sigmaXXMatrix.dat");

        HEAD std::cout << "Long. cond.:" << std::endl;
        HEAD std::cout << "============" << std::endl;
        for(int alpha = 1; alpha <= m_nLayers; alpha++) {
            sigmaXX[alpha] = 0.0;
            for(int beta = 1; beta <= m_nLayers; beta++)
	      if(beta >= alpha) {
                    sigmaXX[alpha] += piXX[alpha][beta];
		    HEAD sigmaMatrix << alpha << " " << beta << " " << piXX[alpha][beta] << std::endl; 
	      }
	      else {
                    sigmaXX[alpha] += piXX[beta][alpha];
                    HEAD sigmaMatrix << alpha << " " << beta << " " << piXX[beta][alpha] << std::endl;
	      }
	    HEAD std::cout << alpha << " " << sigmaXX[alpha] << std::endl;
	    HEAD sigmaFile << alpha << " " << sigmaXX[alpha] << std::endl;
        }
        HEAD std::cout << std::endl;
	
	HEAD sigmaFile.close();
	HEAD sigmaMatrix.close();
    }

    if(m_calculateRhoXY) {
        double sigmaXY[m_nLayers+2];
        double piXY[m_nLayers+2][m_nLayers+2][m_nLayers+2];

        int counter = 1;
        for(int delta = 1; delta <= m_nLayers; delta++)
            for(int alpha = 1; alpha <= m_nLayers; alpha++)
                for(int gamma = delta; gamma <= m_nLayers; gamma++) {
                    piXY[alpha][gamma][delta] = calculatePiXY(alpha, gamma, delta, mpi);
                    HEAD std::cout << alpha << " " << gamma << " " << delta << " " << piXY[alpha][gamma][delta] << " " << counter << "/" << m_nLayers * m_nLayers * (m_nLayers+1) / 2 << std::endl;
                    counter++;
                }
        HEAD std::cout << "Hall cond.:" << std::endl;
        HEAD std::cout << "===========" << std::endl;
        for(int delta = 1; delta <= m_nLayers; delta++) {
            sigmaXY[delta] = 0.0;
            for(int alpha = 1; alpha <= m_nLayers; alpha++)
                for(int gamma = 1; gamma <= m_nLayers; gamma++)
                    if(gamma >= delta)
                        sigmaXY[delta] += piXY[alpha][gamma][delta];
                    else
                        sigmaXY[delta] += piXY[alpha][delta][gamma];
            HEAD std::cout << delta << " " << sigmaXY[delta] << std::endl;
        }
    }

}

void NanoStructure::solve() {
    std::cout << "ECR: " << (m_doECR ? "true" : "false") << std::endl;
    m_iteration = 1;

    if(m_doECR)
        initializePotentials();
    return;


    if(m_doECR) {
        std::string rhoFilename = "rho_" + utils::toString(m_iteration) + ".dat";
        std::ofstream rhoFile(rhoFilename.c_str());
        for(int n = 1; n <= m_nLayers; n++) {
            rhoFile << n << " " << m_rho[n] - m_rhoBackground[n] << std::endl;
        }
        rhoFile.close();

        // updatePotentials
        updatePotentials();
        std::string potFilename = "V_" + utils::toString(m_iteration) + ".dat";
        std::ofstream potFile(potFilename.c_str());
        for(int n = 1; n <= m_nLayers; n++) {
            potFile << n << " " << m_V[n] << std::endl;
        }
        potFile.close();
    }

    while( (m_iteration <= m_maxIterations) && !isConverged() ) {
        std::cout << "Iteration: " << m_iteration << std::endl;
        // calculate local Greensfunctions and occupations
        for(int n = 1; n <= m_nLayers; n++) {
            std::cout << "Calculating local Green's function for layer " << n << std::endl;
            math::CFunction F;
            math::CFunction G;

            // calculat F(w) and G(w)
            F_local(n, F);
            G_local(n, G);

            //update similarity record
            if(m_iteration >= 2)
                m_similarities[n] = math::similarity(G, m_greensFunctions[n]);
            m_greensFunctions[n] = G;
            m_fFunctions[n] = F;
            m_greensFunctions[n].write("g_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");
            m_fFunctions[n].write("f_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");
            if(m_doECR)
                m_rho[n] = calculateOccupationFromG(n);
        }
        for(int n = 1; n <= m_nLayers; n++)
            std::cout << "layer " << n << ": similarity: " << m_similarities[n] << std::endl;

        //calculate effective media
        for(int n = 1; n <= m_nLayers; n++) {
            m_effectiveMedia[n].resize(m_greensFunctions[n].getSize());
            for(int nOmega = 0; nOmega < m_greensFunctions[m_nLayers].getSize(); nOmega++) {
                double omega = m_greensFunctions[n].getArgument(nOmega);
                m_effectiveMedia[n].set(nOmega, omega, -(omega + std::complex<double>(0.0, 1.0) * m_delta) - m_U[n]/2.0 - m_mu[n] + m_V[n]
                                        + m_fFunctions[n].getValue(nOmega) / m_greensFunctions[n].getValue(nOmega) );
            }
            m_effectiveMedia[n].rescale(1.0 / m_scaling);
            m_effectiveMedia[n].write("d_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");
        }

        // do NRG now
        for(int n = 1; n <= m_nLayers; n++) {
            // feed that to the nrg as a hybridization function
            // solve impurity problem
            // update self energy
            config::Configuration& configuration = config::Configuration::getInstance();
            configuration.readFile("config.txt");

            nrg::broadening::FFTBroadener broadener;
            nrg::chain::HybridizationProvider chainProvider;
            chainProvider.setHybridization(m_effectiveMedia[n], m_effectiveMedia[n]);
            chainProvider.setPHsymmetric(true);
            chainProvider.setSZsymmetric(true);

            nrg::NRG nrg(chainProvider, broadener);
            nrg.configure(configuration);
            nrg.setU(m_U[n]/m_scaling);
            nrg.setEpsF(-m_U[n]/(2.0 * m_scaling));
            nrg.setTemperature(m_T);
            nrg.init();
            nrg.showInfo();
            nrg.solve();

            math::CFunction selfEnergy;
            nrg.getSelfEnergy(selfEnergy);
            selfEnergy.rescale(m_scaling);
            m_selfEnergies[n] = selfEnergy;

            m_selfEnergies[n].write("s_" + utils::toString(n) + "_" + utils::toString(m_iteration) + ".dat");


        }
        m_iteration++;
    }


    double sigmaXX[m_nLayers+2];
    double sigmaXY[m_nLayers+2];
    double piXY[m_nLayers+2][m_nLayers+2][m_nLayers+2];

    // post processing

    for(int alpha = 1; alpha <= m_nLayers; alpha++) {
        sigmaXX[alpha] = 0.0;
        for(int beta = 1; beta <= m_nLayers; beta++) {
            double sXX = calculateSigmaXX(alpha, beta);
            std::cout << alpha << " " << beta << " " << sXX << std::endl;
            sigmaXX[alpha] += sXX;
        }
    }
    std::cout << "Conductivities:" << std::endl;
    std::cout << "===============" << std::endl;
    for(int alpha = 1; alpha <= m_nLayers; alpha++)
        std::cout << alpha << " " << sigmaXX[alpha] << std::endl;


    std::cout << "piXy:" << std::endl;
    std::cout << "==================" << std::endl;


    for(int delta = 1; delta <= m_nLayers; delta++) {
        sigmaXY[delta] = 0.0;
        for(int alpha = 1; alpha <= m_nLayers; alpha++) {
            for(int gamma = 1; gamma <= m_nLayers; gamma++) {
                piXY[alpha][gamma][delta] = calculatePiXY(alpha, gamma, delta);
                sigmaXY[delta] += piXY[alpha][gamma][delta];
                std::cout << alpha << " " << gamma << " " << delta << " " << piXY[alpha][gamma][delta] << std::endl;
            }
        }
    }
    std::cout << "piXy:" << std::endl;
    std::cout << "==================" << std::endl;
    for(int delta = 1; delta <= m_nLayers; delta++) {
        std::cout << delta << " " << sigmaXY[delta] << std::endl;
    }

}

void NanoStructure::writeBinary(std::string directory) {
    using namespace utils::io;
    if((directory.size() > 0) && (directory[directory.size()-1] != '/'))
        directory.append(std::string("/"));
    utils::io::createDirectoryStructure(directory);

    std::string filename = directory + "structure.dat";
    std::ofstream outFile(filename.c_str(), std::ios::binary);
    if(!outFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    //global
    binary_write(&outFile, m_nLayers);
    binary_write(&outFile, m_iteration);

    //layers
    for(int n = 0; n < m_nLayers+2; n++) {
        m_selfEnergies[n].writeBinary(directory + utils::toString(n) + ".dat");
        binary_write(&outFile, m_U[n]);
        binary_write(&outFile, m_tL[n]);
        binary_write(&outFile, m_tR[n]);
        binary_write(&outFile, m_rho[n]);
        binary_write(&outFile, m_rhoBackground[n]);
        binary_write(&outFile, m_V[n]);
        binary_write(&outFile, m_mu[n]);
    }
    outFile.close();
}

void NanoStructure::readBinary(std::string directory) {
    using namespace utils::io;
    if((directory.size() > 0) && (directory[directory.size()-1] != '/'))
        directory.append(std::string("/"));

    free();

    std::string filename = directory + "structure.dat";
    std::ifstream inFile(filename.c_str(), std::ios::binary);
    if(!inFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }

    binary_read(&inFile, m_nLayers);
    allocate(m_nLayers);
    binary_read(&inFile, m_iteration);

    for(int n = 0; n < m_nLayers+2; n++) {
        m_selfEnergies[n].readBinary(directory + utils::toString(n) + ".dat");
        binary_read(&inFile, m_U[n]);
        binary_read(&inFile, m_tL[n]);
        binary_read(&inFile, m_tR[n]);
        binary_read(&inFile, m_rho[n]);
        binary_read(&inFile, m_rhoBackground[n]);
        binary_read(&inFile, m_V[n]);
        binary_read(&inFile, m_mu[n]);

    }
    inFile.close();
}

void NanoStructure::initializePotentials() {
    HEAD std::cout << "Initializing potentials" << std::endl;
    for(int iter = 1; iter < 10000; iter++) {
        m_rhoBar = 0.0;
        for(int i = 1; i <= m_nLayers; i++) {
            m_rho[i] = m_n[i].interpolate(m_mu[i] - m_V[i]);
            m_rhoBar += m_rho[i] - m_rhoBackground[i];
        }
        m_rhoBar /= (double) m_nLayers;

        // calculate the potential updates
        for(int i = 1; i <= m_nLayers; i++) {
            double pot = 0.0;
            for(int j = 1; j <= m_nLayers; j++) {
                double distance = std::abs((double) (i-j));
                pot +=  -(m_rho[j] - m_rhoBackground[j] - m_rhoBar) * 0.4 * distance;
            }
            m_V[i] = 0.01 * pot + 0.99 * m_V[i];
        }

    }
    HEAD std::cout << "rhoBar: " << m_rhoBar << std::endl;

    HEAD {
        std::ofstream potFile(std::string("Vi.dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            potFile << i << " " << m_V[i] << " " << std::endl;
        potFile.close();
        std::ofstream occFile(std::string("rhoi.dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            occFile << i << " " << m_rho[i] << " " << m_rho[i]-m_rhoBackground[i] << std::endl;
        occFile.close();
    }
}

void NanoStructure::updatePotentials(mpi::OpenMPI& mpi) {
    // see if we have excess charge on board
    // and lets try to get rid of it by adjusting mu

    //calculate rhoBar
    double rhoBar = 0.0;
    for(int i = 1; i <= m_nLayers; i++) {
        rhoBar += m_rho[i] - m_rhoBackground[i];
    }
    rhoBar /= (double) m_nLayers;

    // calculate the potential updates
    // for the first 20 iterations we let the heterostructure
    // "settle" into the initial potentials
    double potential[m_nLayers+2];
    for(int i = 1; i <= m_nLayers; i++) {
        potential[i] = 0.0;
        for(int j = 1; j <= m_nLayers; j++) {
            double distance = std::abs((double) (i-j));
            potential[i] +=  -(m_rho[j] - m_rhoBackground[j] - rhoBar) * 0.4 * distance;
        }
        if(m_iteration > 20)
            m_V[i] = m_alphaV * potential[i] + (1.0 - m_alphaV) * m_V[i];
    }

    HEAD {
        std::ofstream occFile(std::string("rho_" + utils::toString(m_iteration) + ".dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            occFile << i << " " << m_rho[i]-m_rhoBackground[i] << " " << m_rho[i] << std::endl;
        occFile.close();

        std::ofstream adjustmentFile(std::string("adjustment.dat").c_str(), std::ios::app);
        adjustmentFile << m_iteration << " " << rhoBar <<  std::endl;
        adjustmentFile.close();
        std::ofstream potFile(std::string("V_" + utils::toString(m_iteration) + ".dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            potFile << i << " " << m_V[i] << " " << potential[i] << std::endl;
        potFile.close();
    }

}


void NanoStructure::updatePotentials() {

    //step 0 -- see if we have excess charge on board
    //          and lets try to get rid of it by adjusting mu
    double mu = adjustMu();

    //did this work?



    /*
    // calculate rhoBar
    m_rhoBar = 0.0;
    for(int i = 1; i <= m_nLayers; i++)
        m_rhoBar += m_rho[i] - m_rhoBackground[i];
    m_rhoBar /= ((double) m_nLayers);
    //calculate totalChange
    double totalChange  = 0.0;
    for(int i = 1; i <= m_nLayers; i++) {
        double potential = 0.0;
        for(int j = 1; j <= m_nLayers; j++) {
            double distance = std::abs((int) (i-j));
            potential +=  -(m_rho[j] - m_rhoBackground[j] - m_rhoBar * m_rhoBarFactor ) * distance;
        }
        double oldV = m_V[i];
        totalChange += std::abs((double) (oldV - potential));
    }

    totalChange /= (double)  m_nLayers;

    HEAD    std::cout << "rhoBar: " << m_rhoBar << std::endl;
    HEAD    std::cout << "total change: " << totalChange << std::endl;

    // there is the possibility here to implment some function m_alphaV( totalChange )
    // here to speed up convergence. Dragons ahead though. I am speaking from
    // experience.

    // actually do the update
    double potential[m_nLayers+2];
    for(int i = 1; i <= m_nLayers; i++) {
        potential[i] = 0.0;
        for(int j = 1; j <= m_nLayers; j++) {
            double distance = std::abs((int) (i-j));
            potential[i] +=  -(m_rho[j] - m_rhoBackground[j] - m_rhoBar * m_rhoBarFactor ) * distance;
        }
        double oldV = m_V[i];
        m_V[i] = m_alphaV * potential[i] + (1.0 - m_alphaV) * oldV;
    }
    HEAD {
        std::ofstream rhobarFile(std::string("rhoBar.dat").c_str(), std::ios::app);
        rhobarFile << m_iteration << " " << m_rhoBar << std::endl;
        rhobarFile.close();
        std::ofstream potFile(std::string("V_" + utils::toString(m_iteration) + ".dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            potFile << i << " " << m_V[i] << " " << potential[i] << std::endl;
        potFile.close();
        std::ofstream occFile(std::string("rho_" + utils::toString(m_iteration) + ".dat").c_str());
        for(int i = 1; i <= m_nLayers; i++)
            occFile << i << " " << m_rho[i]-m_rhoBackground[i] << " " << m_rho[i] << std::endl;
        occFile.close();
    }
    */
}

double NanoStructure::adjustMu() {
    double mu_min = -20.0;
    double mu_max = 20.0;
    double mu;
    while(true) {
        mu = (mu_min + mu_max) / 2.0;
        double excess = totalExcessCharge(mu);
        //std::cout << mu << " " << excess << std::endl;
        if(excess < 0)
            mu_min = mu;
        else
            mu_max = mu;
        if((mu_max - mu_min) < 1e-5)
            break;
    }
    return mu;
}

double NanoStructure::totalExcessCharge(double mu) {
    double totalExcess = 0.0;
    for(int n = 1; n <= m_nLayers; n++) {
        double occ = calculateOccupationFromG(n, mu);
        totalExcess += occ - m_rhoBackground[n];
    }
    return totalExcess;
}


double NanoStructure::calculateOccupationFromG(int layer, double mu) {
    assert(layer >= 1);
    assert(layer <= m_nLayers);
    int N = m_greensFunctions[layer].getSize();
    double beta = 1.0 / m_T;
    double integral = 0.0;
    double occIntegral = 0.0;
    for(int n = 0; n < N-1; n++) {
        double avgValue = (m_greensFunctions[layer].getValueImag(n+1) + m_greensFunctions[layer].getValueImag(n))/2.0;
        double avgArg = (m_greensFunctions[layer].getArgument(n+1) + m_greensFunctions[layer].getArgument(n))/2.0;
        double step = m_greensFunctions[layer].getArgument(n+1) - m_greensFunctions[layer].getArgument(n);
        integral += step * avgValue;
        occIntegral += step * avgValue * math::physics::fermi(beta, avgArg-mu);
    }
    double norm = -m_greensFunctions[layer].integrate(m_greensFunctions[layer].getArgument(0), m_greensFunctions[layer].getArgument(N-1)).imag()/M_PI;
    return 2.0 * occIntegral/integral;
}

double NanoStructure::calculateSigmaXX(int alpha, int beta, mpi::OpenMPI& mpi) {
  /*    IntegrationParametersSigmaXX parameters;
    parameters.alpha = alpha;
    parameters.beta = beta;
    parameters.structure = this;
    math::Integral sigmaXXInt;
    sigmaXXInt.setAccuracy(1e-5,1e-5);
    sigmaXXInt.setBounds(3, -4.0 * m_T, 0.0, 4.0 * m_T );
    sigmaXXInt.setIntegralKernels(&NanoStructure::sigmaXX, NULL);
    sigmaXXInt.setIntegrationParameters(&parameters);
    sigmaXXInt.setParameterFunction(m_selfEnergies[alpha]);
    math::CFunction sigmaXXfunc;
    sigmaXXInt.setIntegrals(true, false);
    sigmaXXInt(sigmaXXfunc, mpi);

    // integrate over sigmaXXfunc
    double sigma = 0;
    for(int n = 0; n < sigmaXXfunc.getSize()-1; n++) {
        //std::cout << (sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0 << " " << dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) << std::endl;
        if(std::isnan(dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) * (sigmaXXfunc.getValueReal(n+1) + sigmaXXfunc.getValueReal(n))/2.0 * (sigmaXXfunc.getArgument(n+1) - sigmaXXfunc.getArgument(n))))
            std::cout << n << " " << (sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0 << std::endl;
        sigma += dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) * (sigmaXXfunc.getValueReal(n+1) + sigmaXXfunc.getValueReal(n))/2.0 * (sigmaXXfunc.getArgument(n+1) - sigmaXXfunc.getArgument(n));
    }
    return -sigma;*/
  double sigma = 0;
  IntegrationParametersSigmaXX parameters;
  parameters.alpha = alpha;
  parameters.beta = beta;
  parameters.structure = this;
  math::Integral sigmaXXInt;
  sigmaXXInt.setAccuracy(1e-5,1e-5);
  sigmaXXInt.setBounds(2, -4.0, 4.0);
  sigmaXXInt.setIntegralKernels(&NanoStructure::sigmaXX, NULL);
  sigmaXXInt.setIntegrationParameters(&parameters);
  sigmaXXInt.setParameterFunction(m_selfEnergies[alpha]);
  math::CFunction sigmaXXfunc;
  sigmaXXInt.setIntegrals(true, false);
  sigmaXXInt(sigmaXXfunc, mpi);

  for(int n = 0; n < sigmaXXfunc.getSize()-1; n++) {
    sigma += math::physics::dfermi( 1.0 / getT(), (sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) * (sigmaXXfunc.getValueReal(n+1) + sigmaXXfunc.getValueReal(n))/2.0 * (sigmaXXfunc.getArgument(n+1) - sigmaXXfunc.getArgument(n));
  }
  return -sigma;
}

double NanoStructure::calculateSigmaXX(int alpha, int beta) {
    IntegrationParametersSigmaXX parameters;
    parameters.alpha = alpha;
    parameters.beta = beta;
    parameters.structure = this;
    math::Integral sigmaXXInt;
    sigmaXXInt.setAccuracy(1e-5,1e-5);
    sigmaXXInt.setBounds(3, -4.0 * m_T, 0.0, 4.0 * m_T );
    sigmaXXInt.setIntegralKernels(&NanoStructure::sigmaXX, NULL);
    sigmaXXInt.setIntegrationParameters(&parameters);
    sigmaXXInt.setParameterFunction(m_selfEnergies[alpha]);
    math::CFunction sigmaXXfunc;
    sigmaXXInt.setIntegrals(true, false);
    sigmaXXInt(sigmaXXfunc);

    // integrate over sigmaXXfunc
    double sigma = 0;
    for(int n = 0; n < sigmaXXfunc.getSize()-1; n++) {
        //std::cout << (sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0 << " " << dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) << std::endl;
        if(std::isnan(dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) * (sigmaXXfunc.getValueReal(n+1) + sigmaXXfunc.getValueReal(n))/2.0 * (sigmaXXfunc.getArgument(n+1) - sigmaXXfunc.getArgument(n))))
            std::cout << n << " " << (sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0 << std::endl;
        sigma += dmft::dos::rhoXX((sigmaXXfunc.getArgument(n+1) + sigmaXXfunc.getArgument(n))/2.0) * (sigmaXXfunc.getValueReal(n+1) + sigmaXXfunc.getValueReal(n))/2.0 * (sigmaXXfunc.getArgument(n+1) - sigmaXXfunc.getArgument(n));
    }
    return -sigma;
}

double NanoStructure::calculatePiXY(int alpha, int gamma, int delta, mpi::OpenMPI& mpi) {

    IntegrationParametersPiXY parameters;
    parameters.alpha = alpha;
    parameters.gamma = gamma;
    parameters.delta = delta;
    parameters.structure = this;
    math::Integral piXYInt;
    piXYInt.setAccuracy(1e-5,1e-5);
    piXYInt.setBounds(3, -4.0, 0.0, 4.0 );
    piXYInt.setIntegralKernels(&NanoStructure::piXY, NULL);
    piXYInt.setIntegrationParameters(&parameters);
    piXYInt.setParameterFunction(m_selfEnergies[alpha]);
    math::CFunction piXYfunc;
    piXYInt.setIntegrals(true, false);
    piXYInt(piXYfunc, mpi);

    // integrate over piXYfunc
    double pi = 0;
    for(int n = 0; n < piXYfunc.getSize()-1; n++)
        pi += dmft::dos::rhoXY((piXYfunc.getArgument(n+1) + piXYfunc.getArgument(n))/2.0) * (piXYfunc.getValueReal(n+1) + piXYfunc.getValueReal(n))/2.0 * (piXYfunc.getArgument(n+1) - piXYfunc.getArgument(n));
    return pi;
}

double NanoStructure::calculatePiXY(int alpha, int gamma, int delta) {

    IntegrationParametersPiXY parameters;
    parameters.alpha = alpha;
    parameters.gamma = gamma;
    parameters.delta = delta;
    parameters.structure = this;
    math::Integral piXYInt;
    piXYInt.setAccuracy(1e-5,1e-5);
    piXYInt.setBounds(3, -4.0, 0.0, 4.0 );
    piXYInt.setIntegralKernels(&NanoStructure::piXY, NULL);
    piXYInt.setIntegrationParameters(&parameters);
    piXYInt.setParameterFunction(m_selfEnergies[alpha]);
    math::CFunction piXYfunc;
    piXYInt.setIntegrals(true, false);
    piXYInt(piXYfunc);

    // integrate over piXYfunc
    double pi = 0;
    for(int n = 0; n < piXYfunc.getSize()-1; n++)
        pi += dmft::dos::rhoXY((piXYfunc.getArgument(n+1) + piXYfunc.getArgument(n))/2.0) * (piXYfunc.getValueReal(n+1) + piXYfunc.getValueReal(n))/2.0 * (piXYfunc.getArgument(n+1) - piXYfunc.getArgument(n));
    return pi;
}

bool NanoStructure::areAllSelfEnergySizesEqual() {
    int size = m_selfEnergies[0].getSize();
    for(int i = 1; i < m_nLayers+2; i++)
        if(m_selfEnergies[i].getSize() != size)
            return false;
    return true;
}

bool NanoStructure::isScalingCompatibleWithSelfEnergies() {
    for(int i = 0; i < m_nLayers+2; i++)
        if( (-m_selfEnergies[i].getArgument(0) != m_scaling) || (m_selfEnergies[i].getArgument(m_selfEnergies[i].getSize()-1) != m_scaling) )
            return false;
    return true;
}

void NanoStructure::showInfo() {
    std::cout << "====================" << std::endl;
    std::cout << "** NanoStructure: **" << std::endl;
    std::cout << "====================" << std::endl;
    std::cout << std::endl;
    std::cout << " * symmetric Nanostructure: " << (m_symmetricNano ? "true" : "false") << std::endl;
    std::cout << " * T: " << m_T << std::endl;
    std::cout << " * layers: " << m_nLayers << std::endl;
    std::cout << std::endl;
    std::cout << " - left lead:" << std::endl;
    std::cout << " * U = " << m_U[0] << std::endl;
    std::cout << " * mu = " << m_mu[0] << std::endl;
    std::cout << " * rho = " << m_rho[0] << std::endl;
    std::cout << " * rhoB = " << m_rhoBackground[0] << std::endl;
    std::cout << " * V = " << m_V[0] << std::endl;
    std::cout << " * similarity = " << m_similarities[0] << std::endl;
    std::cout << std::endl;
    for(int i = 1; i <= m_nLayers; i++) {
        std::cout << " - layer " << i << ":" << std::endl;
        std::cout << " * U = " << m_U[i] << std::endl;
        std::cout << " * mu = " << m_mu[i] << std::endl;
        std::cout << " * rho = " << m_rho[i] << std::endl;
        std::cout << " * rhoB = " << m_rhoBackground[i] << std::endl;
        std::cout << " * V = " << m_V[i] << std::endl;
        std::cout << " * similarity = " << m_similarities[i] << std::endl;
        std::cout << std::endl;
    }
    std::cout << " - right lead:" << std::endl;
    std::cout << " * U = " << m_U[m_nLayers+1] << std::endl;
    std::cout << " * mu = " << m_mu[m_nLayers+1] << std::endl;
    std::cout << " * rho = " << m_rho[m_nLayers+1] << std::endl;
    std::cout << " * rhoB = " << m_rhoBackground[m_nLayers+1] << std::endl;
    std::cout << " * V = " << m_V[m_nLayers+1] << std::endl;
    std::cout << " * similarity = " << m_similarities[m_nLayers+1] << std::endl;
    std::cout << std::endl;

}

}
