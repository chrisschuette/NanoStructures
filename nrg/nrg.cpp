#include "nrg.h"
#include "hilbertsubspace.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>

#include "../math/lapack.h"

//#define DEBUG

namespace nrg {
NRG::NRG(ChainProvider& chainProvider, Broadener& broadener)
    : m_chainProvider(chainProvider)
    , m_broadener(broadener)
    , m_G_Up(0)
    , m_G_Down(0)
    , m_F_Up(0)
    , m_F_Down(0)
    , m_epsF(0)
    , m_U(0)
    , m_temperature(1e-20)
    , m_maxIterations(0) //let temperature decide
    , m_nFirstTruncated(-1)
    , m_clusterEnergy(1e-5)
    , m_energyCutOff(8.0)
    , m_maxHSdimension(1000)
{
    dQ[0] = 1;
    dQ[1] = 0;
    dQ[2] = 0;
    dQ[3] = -1;

    dSz[0] = 0;
    dSz[1] = -1;
    dSz[2] = 1;
    dSz[3] = 0;

    signLME[0] = 1;
    signLME[1] = -1;
    signLME[2] = -1;
    signLME[3] = 1;
}

void NRG::configure(config::Configuration& configuration) {
    try {
        m_epsF = configuration.getDouble("NRG.epsF");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_U = configuration.getDouble("NRG.U");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_temperature = configuration.getDouble("NRG.temperature");
    } catch( libconfig::SettingNotFoundException& e ) {}
    try {
        m_clusterEnergy = configuration.getDouble("NRG.clusterEnergy");
    } catch(libconfig::SettingNotFoundException& e) {}
    try {
        m_energyCutOff = configuration.getDouble("NRG.energyCutOff");
    } catch(libconfig::SettingNotFoundException& e) {}
    try {
        m_maxHSdimension = configuration.getInteger("NRG.maxHSdimension");
    } catch(libconfig::SettingNotFoundException& e) {}
}

void NRG::init() {
    if(m_maxIterations == 0) {
        std::cout << " ** Setting max. iteration by temperature. ** " << std::endl;
        m_maxIterations = (int) std::log(m_temperature)/std::log(m_chainProvider.getLambda())*(-2.0) + 15;
    }

    if( (std::log(m_temperature)/std::log(m_chainProvider.getLambda())*(-2.0) + 9) >= m_maxIterations)
        std::cout << "You are not doing enough NRG iterations." << std::endl;

    // consistency check
    if((std::abs(m_U + 2.0 * m_epsF) > 1e-20) && m_chainProvider.isPHsymmetric()) {
        std::cerr << "Problem not PH symmetric but have PH symmetry enforced." << std::endl;
        throw std::exception();
    }

    //always do an even number of iterations
    m_maxIterations += m_maxIterations%2;

    m_chainProvider.setLength(m_maxIterations);
    m_chainProvider.buildChain();
    //m_chainProvider.dump();

    // initialize the broadenner
    m_broadener.setTemperature(m_temperature);

    // initialize the broadeners for the various correlation functions
    m_G_Up = m_broadener.clone();
    m_G_Up->init();
    m_G_Down = m_broadener.clone();
    m_G_Down->init();

    m_F_Up = m_broadener.clone();
    m_F_Up->init();
    m_F_Down = m_broadener.clone();
    m_F_Down->init();

    setupInitialState();
}


void NRG::setupInitialState() {
    double pref = 2.0 / (1.0 + m_chainProvider.getLambda());
    double e_d = pref * m_epsF;
    double U = pref * m_U;

    HilbertSubSpace* H1 = new HilbertSubSpace(-1, 0); // no electron
    HilbertSubSpace* H2 = new HilbertSubSpace(0, 1);  // spin up electron
    HilbertSubSpace* H3 = new HilbertSubSpace(0, -1); // spin down electron
    HilbertSubSpace* H4 = new HilbertSubSpace(1, 0);  // two electrons

    // make space for eigenenergies
    H1->getEnergies().resize(1);
    H2->getEnergies().resize(1);
    H3->getEnergies().resize(1);
    H4->getEnergies().resize(1);

    // set the initial Eigenenergies
    H1->getEnergies().set(0,0.0);
    H2->getEnergies().set(0,e_d);
    H3->getEnergies().set(0,e_d);
    H4->getEnergies().set(0,2.0 * e_d + U);

    // in the -1 iteration we keep all the states
    H1->setKeptStates(1);
    H2->setKeptStates(1);
    H3->setKeptStates(1);
    H4->setKeptStates(1);

    // set initial matrix elements
    {
        // \hat{c}_{N+1,\uparrow}
        H2->getChainElementsUp().resize(1,1);
        H2->getChainElementsUp().set(0,0,1);
        H4->getChainElementsUp().resize(1,1);
        H4->getChainElementsUp().set(0,0,1);

        H2->getLocalMatrixElementUp().resize(1,1);
        H2->getLocalMatrixElementUp().set(0,0,1);
        H4->getLocalMatrixElementUp().resize(1,1);
        H4->getLocalMatrixElementUp().set(0,0,1);

        H2->getLocalMatrixElementUp2().resize(1,1);
        H2->getLocalMatrixElementUp2().set(0,0,0);
        H4->getLocalMatrixElementUp2().resize(1,1);
        H4->getLocalMatrixElementUp2().set(0,0,1);


        // \hat{c}_{N+1,\downarrow}
        H3->getChainElementsDown().resize(1,1);
        H3->getChainElementsDown().set(0,0,1);
        H4->getChainElementsDown().resize(1,1);
        H4->getChainElementsDown().set(0,0,-1);

        H3->getLocalMatrixElementDown().resize(1,1);
        H3->getLocalMatrixElementDown().set(0,0,1);
        H4->getLocalMatrixElementDown().resize(1,1);
        H4->getLocalMatrixElementDown().set(0,0,-1);

        H3->getLocalMatrixElementDown2().resize(1,1);
        H3->getLocalMatrixElementDown2().set(0,0,0);
        H4->getLocalMatrixElementDown2().resize(1,1);
        H4->getLocalMatrixElementDown2().set(0,0,-1);
    }

    // add them to the HilbertSpace table
    m_hilbertSpaces.setHS(-1, -1, 0, H1);
    m_hilbertSpaces.setHS(-1, 0, 1, H2);
    m_hilbertSpaces.setHS(-1, 0, -1, H3);
    m_hilbertSpaces.setHS(-1, 1, 0, H4);
}

void NRG::setupHamiltonian(int n) {
    double lambda = m_chainProvider.getLambda();

    // get hopping elements and on-site energies from the chain
    double eUp =  m_chainProvider.getOrbitalUp(n).second;
    double tUp =  m_chainProvider.getOrbitalUp(n).first;
    double eDown =  m_chainProvider.getOrbitalDown(n).second;
    double tDown =  m_chainProvider.getOrbitalDown(n).first;

    // maximum possible q value
    // Q is measured relative to half-filling
    int Qmax = n + 2;

    // clean out old energies
    m_energies.clear();

    // go over all HilbertSubSpaces
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            //find the relevant HilbertSubSpaces of the last iteration;
            HilbertSubSpace* H_last[4];
            H_last[0] = m_hilbertSpaces.getHS(n-1, Q+1, Sz);
            H_last[1] = m_hilbertSpaces.getHS(n-1, Q, Sz - 1);
            H_last[2] = m_hilbertSpaces.getHS(n-1, Q, Sz + 1);
            H_last[3] = m_hilbertSpaces.getHS(n-1, Q-1, Sz);

            int dim = 0;
            int keptStates[4];

            for(int k = 0; k < 4; k++) {
                if(H_last[k] != NULL)
                    keptStates[k] = H_last[k]->getKeptStates();
                else
                    keptStates[k] = 0;
                dim += keptStates[k];
            }

            if(dim > 0)
            {
                HilbertSubSpace* HS = new HilbertSubSpace(Q,Sz);
                math::Matrix& H = HS->getHamiltonian();
                HS->setKeptStates(dim);

                for(int k = 0; k < 4; k++)
                    HS->setR(k,keptStates[k]);

                //reserve memory for the Hamiltonian
                H.resize(dim, dim);
                H.zero();

                // build the Hamiltonian
                if((keptStates[0] != 0) && (keptStates[1] != 0))
                    H.copy(H_last[0]->getChainElementsUp(), 0, keptStates[0], keptStates[0], keptStates[1], tUp);
                if((keptStates[0] != 0) && (keptStates[2] != 0))
                    H.copy(H_last[0]->getChainElementsDown(), 0, (keptStates[0] + keptStates[1]), keptStates[0], keptStates[2], tDown);
                if((keptStates[1] != 0) && (keptStates[3] != 0))
                    H.copy(H_last[1]->getChainElementsDown(), keptStates[0], keptStates[0] + keptStates[1] + keptStates[2], keptStates[1], keptStates[3], tDown);
                if((keptStates[2] != 0) && (keptStates[3] != 0))
                    H.copy(H_last[2]->getChainElementsUp(), keptStates[0] + keptStates[1], keptStates[0] + keptStates[1] + keptStates[2], keptStates[2], keptStates[3], -tUp);

                double _sign_d_Up[4];
                _sign_d_Up[0] = 0.0;
                _sign_d_Up[1] = 1.0;
                _sign_d_Up[2] = 0.0;
                _sign_d_Up[3] = 1.0;

                double _sign_d_Down[4];
                _sign_d_Down[0] = 0.0;
                _sign_d_Down[1] = 0.0;
                _sign_d_Down[2] = 1.0;
                _sign_d_Down[3] = 1.0;

                // place scaled eigenenergies of last iteration on diagonal
                for(int k = 0; k < 4; k++)
                    if(keptStates[k] != 0)
                        for(int i = 0; i < keptStates[k]; i++)
                        {
                            int offset = 0;
                            for(int j = 0; j < k; j++)
                                offset += keptStates[j];
                            H.set(offset+i, offset+i, sqrt(lambda) * H_last[k]->getEnergies().get(i) + _sign_d_Up[k] * eUp + _sign_d_Down[k] * eDown);
                        }

                //diagonalize the Hamiltonian
                char jobz, uplo;
                int lda, lwork, info;
                double *work;

                jobz = 'V';
                uplo = 'U';
                lda = dim;
                lwork = 3*dim-1;
                work = new double[lwork];
                HS->getEnergies().resize(dim);

                lapack::dsyev_(&jobz, &uplo, &dim, HS->getHamiltonian().getData(), &lda, HS->getEnergies().getData(), work, &lwork, &info);
                delete [] work;

                //move the new eigenenergies into the energies array
                unsigned int currentSize = m_energies.size();
                m_energies.resize( currentSize + dim);
                std::copy( HS->getEnergies().getData(), HS->getEnergies().getData() + dim, m_energies.begin() + currentSize );

                // save to HilbertSpaces table
                m_hilbertSpaces.setHS(n, Q, Sz, HS);
            }
        }
    }


    //correct the energies in the subspaces
    std::sort(m_energies.begin(), m_energies.end());
    double gs_energy = m_energies[0];
    for(int i = 0; i < (int) m_energies.size(); i++)
    {
        m_energies[i] = m_energies[i] - gs_energy;
    }

    // go over all HilbertSubSpaces
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q ,Sz);
            if(HS != NULL)
                for(int j = 0; j < HS->getEnergies().getLength(); j++)
                    HS->getEnergies().set(j, HS->getEnergies().get(j) - gs_energy);
        }
    }
    //clean up
    //discard old matrix elements.
    deleteChainOperatorElements(n-1);
    deleteChainOperatorElements(n-1);
}

void NRG::truncateStates(int n) {
    int last_iteration = m_maxIterations - 1;
    int Qmax = n + 2;

    /*
     *  SYMMETRIES: SZ
     */

    if(m_chainProvider.isSZsymmetric())
    {
        //  std::cout << "SYMMETRY: Enforcing SZ symmetry." << std::endl;

        //symmetrize energies regarding sz
        //Q loop
        for(int Q = Qmax; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int Sz = 0; Sz <= Szmax; Sz++)
            {
                // get the two HilberSubSpaces that should be averaged.
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q ,Sz);
                HilbertSubSpace* HS_minusSz = m_hilbertSpaces.getHS(n, Q ,-Sz);
                if((HS != NULL) && (HS_minusSz != NULL)) {
                    if(HS->getEnergies().getLength() == HS_minusSz->getEnergies().getLength())
                    {
                        for(int i = 0; i < HS->getEnergies().getLength(); i++)
                        {
                            double average = (HS->getEnergies().get(i) + HS_minusSz->getEnergies().get(i))/2.0;
                            HS->getEnergies().set(i,average);
                            HS_minusSz->getEnergies().set(i,average);

                        }
                    }
                    else
                    {
                        std::cerr << "WARNING! CAN'T SYMMETRIZE: Q=" << Q << " Sz=" << Sz << std::endl;
                    }
                }
            }
        }
    }

    /*
      SYMMETRIES: PH
      */

    if(m_chainProvider.isPHsymmetric())
    {
        //  std::cout << "SYMMETRY: Enforcing PH symmetry." << std::endl;
        //symmetrize energies regarding ph
        //Q loop
        for(int Q = 0; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int Sz = -Szmax; Sz <= Szmax; Sz++)
            {
                // get the two HilberSubSpaces that should be averaged.
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q ,Sz);
                HilbertSubSpace* HS_minusQ = m_hilbertSpaces.getHS(n, -Q ,Sz);
                if((HS != NULL) && (HS_minusQ != NULL)) {
                    if(HS->getEnergies().getLength() == HS_minusQ->getEnergies().getLength())
                    {
                        for(int i = 0; i < HS->getEnergies().getLength(); i++)
                        {
                            double average = (HS->getEnergies().get(i) + HS_minusQ->getEnergies().get(i))/2.0;
                            HS->getEnergies().set(i,average);
                            HS_minusQ->getEnergies().set(i,average);

                        }
                    }
                    else
                    {
                        std::cerr << "WARNING! CAN'T SYMMETRIZE: Q=" << Q << " Sz=" << Sz << std::endl;
                    }
                }
            }
        }

    }

    if(n == last_iteration)
    {
        //discard all states
        int Qmax = n + 2;

        //Q loop
        for(int Q = Qmax; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int Sz = -Szmax; Sz <= Szmax; Sz++)
            {
                //set kept states
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q ,Sz);
                if(HS != NULL)
                    HS->setKeptStates(0);
            }
        }
        return;
    }

    //find index of first high energy state
    int highEnergyIndex = 0;

    //    std::cout << "number of energies: " << m_energies.size() << std::endl;
    //    std::cout << "max energy: " << m_energies[m_energies.size() - 1] << std::endl;

    // find the first energy > cutoff.
    while((highEnergyIndex < (int) m_energies.size()) && (m_energies[highEnergyIndex] < m_energyCutOff))
        highEnergyIndex++;

    // keep at most m_maxHSDim states
    highEnergyIndex = std::min(m_maxHSdimension, highEnergyIndex);

    if(highEnergyIndex < (int) m_energies.size())
    {
        if(m_nFirstTruncated == -1)
        {
            //   std::cout << "State truncation starts in iteration " << n << "." << std::endl;
            m_nFirstTruncated = n;
        }
        //        std::cout << "First high energy state: " << m_energies[highEnergyIndex] << std::endl;
    }

    // the goal is to set the variable E_threshold to a value
    // between the last energy in the energy cluster and the
    // first energy above that

    int index;
    for(index = highEnergyIndex; index < (int) m_energies.size(); index++) {
        double e_diff = m_energies[index] - m_energies[index-1];
        if(e_diff > m_clusterEnergy)
            break;
    }

    if(index > 1.2 *  m_maxHSdimension) {
        std::cout << "This is getting out of hand - reduce cluster energy." << std::endl;
        throw std::exception();
    }

    double E_threshold;

    if(index >= (int) m_energies.size())
        E_threshold = m_energies[m_energies.size()-1] + 0.1; // nothing gets discarded
    else
        E_threshold = (m_energies[index] + m_energies[index-1])/2.0;

    //    std::cout << "Truncating everything above E_threshold = " << E_threshold << std::endl;

    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q ,Sz);
            if(HS != NULL)
            {
                int i = 0;
                while((i < HS->getEnergies().getLength()) && (HS->getEnergies().get(i) < E_threshold) )
                    i++;

                //set kept states
                HS->setKeptStates(i);
            }
        }
    }
}


void NRG::deleteChainOperatorElements(int iteration) {
    // don't clean iteration -1
    if(iteration < 0)
        return;
    int Qmax = iteration + 2;

    // iterate over HilbertSpace
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q ,Sz);
            if(HS != NULL) {
                HS->getChainElementsUp().resize(0,0);
                HS->getChainElementsDown().resize(0,0);
            }
        }
    }
}

void NRG::deleteTransformationMatrices(int iteration) {
    // don't clean iteration -1
    if(iteration < 0)
        return;
    int Qmax = iteration + 2;

    // iterate over HilbertSpace
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q ,Sz);
            if(HS != NULL) {
                HS->getHamiltonian().resize(0,0);
            }
        }
    }
}

void NRG::deleteDensityMatrices(int iteration) {
    // don't clean iteration -1
    if(iteration < 0)
        return;
    int Qmax = iteration + 2;

    // iterate over HilbertSpace
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q ,Sz);
            if(HS != NULL) {
                HS->getDensityMatrix().resize(0,0);
                HS->getDensityMatrixEigenBasis().resize(0,0);
            }
        }
    }
}

void NRG::deleteImpurityMatrixElements(int iteration) {
    // don't clean iteration -1
    if(iteration < 0)
        return;
    int Qmax = iteration + 2;

    // iterate over HilbertSpace
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int Sz = -Szmax; Sz <= Szmax; Sz++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q ,Sz);
            if(HS != NULL) {
                HS->getLocalMatrixElementDown().resize(0,0);
                HS->getLocalMatrixElementDown2().resize(0,0);
                HS->getLocalMatrixElementUp().resize(0,0);
                HS->getLocalMatrixElementUp2().resize(0,0);            }
        }
    }
}

void NRG::propagateChainOperatorElementsUp(int iteration) {
    int Qmax = iteration + 2;
    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(iteration, Q - 1, SZ - 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getKeptStates();
                    off_target_rows = 0;
                    target_cols = H_other->getKeptStates();
                    off_target_cols = 0;
                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getChainElementsUp().resize(target_rows, target_cols);
                        HS->getChainElementsUp().zero();

                        /*
                         *  i = 1 && j = 0
                         */

                        int offset_i = HS->getR(0);
                        int offset_j = 0;

                        // variables for matrix multiplication
                        double one = 1.0;
                        char N = 'N';
                        char T = 'T';

                        int k_max = HS->getR(1);
                        if(k_max > 0)
                        {
                            //copy U matrizes
                            math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(1), target_rows);
                            math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(0), target_cols);
                            double factor = 1.0;
                            lapack::dgemm_(&T, &N, &target_rows, &target_cols, &k_max, &factor, UT->getData(), &k_max, U->getData(), &k_max, &one, HS->getChainElementsUp().getData(), &target_rows);
                            delete U;
                            delete UT;
                        }

                        /*
                         *  i = 3 && j =2
                         */
                        offset_i = HS->getR(0) + HS->getR(1) + HS->getR(2);
                        offset_j = H_other->getR(0) + H_other->getR(1);

                        k_max = HS->getR(3);
                        if(k_max > 0)
                        {
                            //copy U matrizes
                            math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(3), target_rows);
                            math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(2), target_cols);

                            double factor = 1.0;
                            lapack::dgemm_(&T, &N, &target_rows, &target_cols, &k_max, &factor, UT->getData(), &k_max, U->getData(), &k_max, &one, HS->getChainElementsUp().getData(), &target_rows);
                            delete U;
                            delete UT;
                        }
                    }
                }
            }
        }
    }
}

void NRG::propagateChainOperatorElementsDown(int iteration) {
    int Qmax = iteration + 2;
    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(iteration, Q - 1, SZ + 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getKeptStates();
                    off_target_rows = 0;
                    target_cols = H_other->getKeptStates();
                    off_target_cols = 0;

                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getChainElementsDown().resize(target_rows, target_cols);
                        HS->getChainElementsDown().zero();

                        /*
                         * i = 2 && j = 0
                         */

                        int offset_i = HS->getR(0) + HS->getR(1);
                        int offset_j = 0;

                        // variables for matrix multiplication
                        double one = 1.0;
                        char N = 'N';
                        char T = 'T';

                        int k_max = HS->getR(2);
                        if(k_max > 0)
                        {
                            //copy U matrizes
                            math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(2), target_rows);
                            math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(0), target_cols);

                            double factor = 1.0;
                            lapack::dgemm_(&T, &N, &target_rows, &target_cols, &k_max, &factor, UT->getData(), &k_max, U->getData(), &k_max, &one, HS->getChainElementsDown().getData(), &target_rows);
                            delete U;
                            delete UT;
                        }

                        /*
                         * i = 3 && j = 1
                         */

                        offset_i = HS->getR(0) + HS->getR(1) + HS->getR(2);
                        offset_j = H_other->getR(0);

                        k_max = HS->getR(3);
                        if(k_max > 0)
                        {
                            //copy U matrizes
                            math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(3), target_rows);
                            math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(1), target_cols);

                            double factor = -1;
                            lapack::dgemm_(&T, &N, &target_rows, &target_cols, &k_max, &factor, UT->getData(), &k_max, U->getData(), &k_max, &one, HS->getChainElementsDown().getData(), &target_rows);
                            delete U;
                            delete UT;
                        }
                    }
                }
            }
        }
    }
}

void NRG::builDM() {
    int last_iteration = m_maxIterations-1;
    double beta = 1.0/m_temperature;

    double Znd[m_maxIterations];
    double Ztot = 0;

    // rescale (discarded) energies to absolute units & calculate partition function of discarded states (iteration-wise)
    for(int iteration = 0; iteration <= last_iteration; iteration++)
    {
        Znd[iteration] = 0.0;
        double factor = 0.5 * (1.0 + 1.0 / m_chainProvider.getLambda()) * std::pow(m_chainProvider.getLambda(), -0.5*((double) (iteration - 1)));
        int Qmax = iteration + 2;
        //Q loop
        for(int Q = Qmax; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int SZ = -Szmax; SZ <= Szmax; SZ++)
            {
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q, SZ);
                if(HS != NULL)
                {
                    //rescale all energies
                    for(int j = 0; j < HS->getEnergies().getLength(); j++)
                        HS->getEnergies().set(j,factor * HS->getEnergies().get(j));
                    for(int j = HS->getKeptStates(); j < HS->getEnergies().getLength(); j++)
                    {
                        double znd_contrib = exp(-beta * HS->getEnergies().get(j));
                        assert(znd_contrib != std::numeric_limits<double>::infinity());
                        assert(!std::isnan(znd_contrib));
                        Znd[iteration] += znd_contrib;
                    }
                }
            }
        }
        double Z = Znd[iteration];
        double Lambda = std::pow(4.0, (last_iteration - iteration));
        assert(!std::isnan(Z));
        assert(!std::isnan(Lambda));
        assert(Z != std::numeric_limits<double>::infinity());
        assert(Lambda != std::numeric_limits<double>::infinity());
        Ztot += Z * Lambda;
    }
    //   std::cout << "Ztot: " << Ztot << std::endl;

    double wn[m_maxIterations];

#ifdef DEBUG
    std::ofstream weightsFile("debug_weights.txt");
#endif
    //calculate weights
    double sum = 0.0;
    for(int iteration = 0; iteration <= last_iteration; iteration++)
    {
        double Z = Znd[iteration];
        double Lambda = std::pow(4.0, (last_iteration - iteration));
        wn[iteration] = Z * Lambda / Ztot;
#ifdef DEBUG
        weightsFile << iteration << " " << wn[iteration] << std::endl;
#endif
        sum += wn[iteration];
    }
    //    std::cout << "sum: " << sum << std::endl;
#ifdef DEBUG
    weightsFile.close();
#endif

    //    std::cout << "Building density matrix...";

    // construct DM for the last_iteration
    int iteration = last_iteration;

    int Qmax = iteration + 2;

    double weight = 1.0 / Ztot;

    double trace = 0.0;

    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q, SZ);
            if(HS != NULL)
            {
                //reserve memory for the Density Matrix (total number of states in last iteration)
                HS->getDensityMatrix().resize(HS->getEnergies().getLength(), HS->getEnergies().getLength());
                HS->getDensityMatrix().zero();

                //reserve memory for the Density Matrix (total number of states in last iteration)
                HS->getDensityMatrixEigenBasis().resize(HS->getEnergies().getLength(), HS->getEnergies().getLength());
                HS->getDensityMatrixEigenBasis().zero();

                for(int i = 0; i < HS->getEnergies().getLength(); i++)
                {
                    HS->getDensityMatrixEigenBasis().set(i,i,weight * exp(-beta * HS->getEnergies().get(i)));
                    trace += HS->getDensityMatrixEigenBasis().get(i,i);
                }

                // project it back!
                math::Matrix temp(HS->getEnergies().getLength(), HS->getEnergies().getLength());

                int n = HS->getEnergies().getLength();
                double zero = 0.0;
                double one = 1.0;
                char N = 'N';
                char T = 'T';
                // rho * U^T
                lapack::dgemm_(&N, &T, &n, &n, &n, &one, HS->getDensityMatrixEigenBasis().getData(), &n, HS->getHamiltonian().getData(), &n, &zero, temp.getData(), &n);
                // U * (rho * U^T)
                lapack::dgemm_(&N, &N, &n, &n, &n, &one, HS->getHamiltonian().getData(), &n, temp.getData(), &n, &zero, HS->getDensityMatrix().getData(), &n);
            }
        }
    }

    int active;

    for(int iteration = last_iteration - 1; iteration >= m_nFirstTruncated; iteration--)
    {
        trace = 0.0;
        weight = std::pow(4.0, last_iteration - iteration) / Ztot;

        int Qmax = iteration + 2;
        for(int Q = Qmax; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int SZ = -Szmax; SZ <= Szmax; SZ++)
            {
                active = 0;
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration, Q, SZ);
                if(HS != NULL)
                {
                    HilbertSubSpace* H_next[4];
                    H_next[0] = m_hilbertSpaces.getHS(iteration+1, Q+1, SZ);
                    H_next[1] = m_hilbertSpaces.getHS(iteration+1, Q, SZ-1);
                    H_next[2] = m_hilbertSpaces.getHS(iteration+1, Q, SZ+1);
                    H_next[3] = m_hilbertSpaces.getHS(iteration+1, Q-1, SZ);

                    int size_kept_part = 0;

                    for(int hi = 0; hi < 4; hi++)
                        if((H_next[hi] != NULL) && (H_next[hi]->getR(3 - hi) > 0))
                            active++;


                    //TODO: is this really correct????
                    if(active == 4)
                    {
                        int getR = H_next[0]->getR(3);
#ifdef DEBUG
                        //check dimensionality
                        for(int i2 = 0; i2 < 4; i2++)
                        {
                            if(H_next[i2]->getR(3 - i2) != getR)
                            {
                                std::cout << "Q:" << Q << " Sz: " << SZ << " iteration: " << iteration << " active = " << active<< std::endl;

                                std::cout << "ERROR!" << std::endl;
                                for(int i2 = 0; i2 < 4; i2++)
                                {
                                    std::cout << i2 << " " << H_next[i2]->getR(3 - i2) << std::endl;
                                }
                                return;
                            }
                        }
#endif
                        if(HS->getKeptStates() != getR)
                            std::cout << "ERROR!" << std::endl;
                        size_kept_part = HS->getKeptStates();
                    }
                    else if(active == 0) {
                        size_kept_part = 0;
                    }
                    else
                    {
                        std::cout << "active != 4 && active != 0" << std::endl;
                        throw std::exception();
                    }

                    int size_discarded_part = HS->getEnergies().getLength() - HS->getKeptStates();

                    HS->getDensityMatrix().resize(size_kept_part + size_discarded_part, size_kept_part + size_discarded_part);
                    //		std::cout << "Reserving " << size_kept_part + size_discarded_part << " x " << size_kept_part + size_discarded_part << std::endl;

                    HS->getDensityMatrix().zero();

                    HS->getDensityMatrixEigenBasis().resize(size_kept_part + size_discarded_part, size_kept_part + size_discarded_part);
                    //		std::cout << "Reserving " << size_kept_part + size_discarded_part << " x " << size_kept_part + size_discarded_part << std::endl;

                    HS->getDensityMatrixEigenBasis().zero();

                    if(size_kept_part > 0)
                    {
                        int offset[4];
                        offset[0] = H_next[0]->getR(0) + H_next[0]->getR(1) + H_next[0]->getR(2);
                        offset[1] = H_next[1]->getR(0) + H_next[1]->getR(1);
                        offset[2] = H_next[2]->getR(0);
                        offset[3] = 0;

                        for(int r1 = 0; r1 < size_kept_part; r1++)
                        {
                            for(int c1 = 0; c1 < size_kept_part; c1++)
                            {
                                double value = H_next[0]->getDensityMatrix().get(offset[0] + r1, offset[0] + c1);
                                value += H_next[1]->getDensityMatrix().get(offset[1] + r1, offset[1] + c1);
                                value += H_next[2]->getDensityMatrix().get(offset[2] + r1, offset[2] + c1);
                                value += H_next[3]->getDensityMatrix().get(offset[3] + r1, offset[3] + c1);
                                HS->getDensityMatrixEigenBasis().set(r1, c1, value);
                            }
                        }
                    }

                    // discarded diagonal part
                    for(int d1 = size_kept_part; d1 < HS->getEnergies().getLength(); d1++)
                    {
                        double value = weight * exp(-beta * HS->getEnergies().get(d1));
                        HS->getDensityMatrixEigenBasis().set(d1,d1,value);
                    }

                    //calculate trace
                    for(int t = 0; t < HS->getEnergies().getLength(); t++)
                        trace += HS->getDensityMatrixEigenBasis().get(t,t);

                    //project it back!
                    math::Matrix temp(HS->getEnergies().getLength(), HS->getEnergies().getLength());

                    int n = HS->getEnergies().getLength();
                    double zero = 0.0;
                    double one = 1.0;
                    char N = 'N';
                    char T = 'T';

                    // rho * U^T
                    lapack::dgemm_(&N, &T, &n, &n, &n, &one, HS->getDensityMatrixEigenBasis().getData(), &n, HS->getHamiltonian().getData(), &n, &zero, temp.getData(), &n);
                    // U * (rho * U^T)
                    lapack::dgemm_(&N, &N, &n, &n, &n, &one, HS->getHamiltonian().getData(), &n, temp.getData(), &n, &zero, HS->getDensityMatrix().getData(), &n);
                }
            }
        }
        //clean old DM's
        Qmax = (iteration+1) + 2;
        for(int Q = Qmax; Q >= -Qmax; Q--)
        {
            int Szmax = (Qmax - abs(Q));
            for(int SZ = -Szmax; SZ <= Szmax; SZ++)
            {
                HilbertSubSpace* HS = m_hilbertSpaces.getHS(iteration+1, Q, SZ);
                if(HS != NULL)
                    HS->getDensityMatrix().resize(0,0);
            }
        }
    }

    //clean old DM's of first Truncated
    Qmax = (m_nFirstTruncated) + 2;
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(m_nFirstTruncated, Q, SZ);
            if(HS != NULL)
                HS->getDensityMatrix().resize(0,0);
        }
    }

    // std::cout << "done." << std::endl;
}

void NRG::propagateLocalMatrixElementDown(int n) {
    int Qmax = n + 2;

    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q - 1, SZ + 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getEnergies().getLength();
                    off_target_rows = 0;
                    target_cols = H_other->getEnergies().getLength();
                    off_target_cols = 0;

                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getLocalMatrixElementDown().resize(target_rows, target_cols);
                        HS->getLocalMatrixElementDown().zero();

                        for(int i = 0; i < 4; i++)
                        {
                            int offset_i = 0;
                            int offset_j = 0;
                            for(int g = 0; g < i; g++) {
                                offset_i += HS->getR(g);
                                offset_j += H_other->getR(g);
                            }

                            // variables for matrix multiplication
                            double one = 1.0;
                            double zero = 0.0;
                            char N = 'N';
                            char T = 'T';

                            //3 Matrizes
                            int k_max = HS->getR(i);
                            int kb_max = H_other->getR(i);
                            HilbertSubSpace* H_last = m_hilbertSpaces.getHS(n - 1, Q + dQ[i], SZ + dSz[i]);
                            if((H_last != NULL) && (k_max > 0) && (kb_max > 0))
                            {
                                //copy U matrizes
                                math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(i), target_rows);
                                math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(i), target_cols);

                                math::Matrix * imp_kept = H_last->getLocalMatrixElementDown().crop(0,0,HS->getR(i), H_other->getR(i));

                                math::Matrix temp;
                                temp.resize(target_rows, kb_max );

                                double factor = signLME[i];

                                lapack::dgemm_(&T, &N, &target_rows, &kb_max, &k_max, &factor, UT->getData(), &k_max, imp_kept->getData() , &k_max, &zero, temp.getData(), &target_rows);
                                lapack::dgemm_(&N, &N, &target_rows, &target_cols, &kb_max, &one, temp.getData(), &target_rows, U->getData(), &kb_max, &one, HS->getLocalMatrixElementDown().getData(), &target_rows);

                                delete imp_kept;

                                delete U;
                                delete UT;
                            }
                        }
                    }
                }
            }
        }
    }
}

void NRG::propagateLocalMatrixElementUp2(int n) {
    int Qmax = n + 2;

    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q - 1, SZ - 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getEnergies().getLength();
                    off_target_rows = 0;
                    target_cols = H_other->getEnergies().getLength();
                    off_target_cols = 0;

                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getLocalMatrixElementUp2().resize(target_rows, target_cols);
                        HS->getLocalMatrixElementUp2().zero();

                        for(int i = 0; i < 4; i++)
                        {
                            int offset_i = 0;
                            int offset_j = 0;
                            for(int g = 0; g < i; g++) {
                                offset_i += HS->getR(g);
                                offset_j += H_other->getR(g);
                            }

                            // variables for matrix multiplication
                            double one = 1.0;
                            double zero = 0.0;
                            char N = 'N';
                            char T = 'T';

                            //3 Matrizes
                            int k_max = HS->getR(i);
                            int kb_max = H_other->getR(i);
                            HilbertSubSpace* H_last = m_hilbertSpaces.getHS(n - 1, Q + dQ[i], SZ + dSz[i]);
                            if((H_last != NULL) && (k_max > 0) && (kb_max > 0))
                            {
                                //copy U matrizes
                                math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(i), target_rows);
                                math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(i), target_cols);

                                math::Matrix * imp_kept = H_last->getLocalMatrixElementUp2().crop(0,0,HS->getR(i), H_other->getR(i));

                                math::Matrix temp;
                                temp.resize(target_rows, kb_max );

                                double factor = signLME[i];

                                lapack::dgemm_(&T, &N, &target_rows, &kb_max, &k_max, &factor, UT->getData(), &k_max, imp_kept->getData() , &k_max, &zero, temp.getData(), &target_rows);
                                lapack::dgemm_(&N, &N, &target_rows, &target_cols, &kb_max, &one, temp.getData(), &target_rows, U->getData(), &kb_max, &one, HS->getLocalMatrixElementUp2().getData(), &target_rows);

                                delete imp_kept;

                                delete U;
                                delete UT;
                            }
                        }
                    }
                }
            }
        }
    }
}

void NRG::propagateLocalMatrixElementDown2(int n) {
    int Qmax = n + 2;

    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q - 1, SZ + 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getEnergies().getLength();
                    off_target_rows = 0;
                    target_cols = H_other->getEnergies().getLength();
                    off_target_cols = 0;

                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getLocalMatrixElementDown2().resize(target_rows, target_cols);
                        HS->getLocalMatrixElementDown2().zero();

                        for(int i = 0; i < 4; i++)
                        {
                            int offset_i = 0;
                            int offset_j = 0;
                            for(int g = 0; g < i; g++) {
                                offset_i += HS->getR(g);
                                offset_j += H_other->getR(g);
                            }

                            // variables for matrix multiplication
                            double one = 1.0;
                            double zero = 0.0;
                            char N = 'N';
                            char T = 'T';

                            //3 Matrizes
                            int k_max = HS->getR(i);
                            int kb_max = H_other->getR(i);
                            HilbertSubSpace* H_last = m_hilbertSpaces.getHS(n - 1, Q + dQ[i], SZ + dSz[i]);
                            if((H_last != NULL) && (k_max > 0) && (kb_max > 0))
                            {
                                //copy U matrizes
                                math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(i), target_rows);
                                math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(i), target_cols);

                                math::Matrix * imp_kept = H_last->getLocalMatrixElementDown2().crop(0,0,HS->getR(i), H_other->getR(i));

                                math::Matrix temp;
                                temp.resize(target_rows, kb_max );

                                double factor = signLME[i];

                                lapack::dgemm_(&T, &N, &target_rows, &kb_max, &k_max, &factor, UT->getData(), &k_max, imp_kept->getData() , &k_max, &zero, temp.getData(), &target_rows);
                                lapack::dgemm_(&N, &N, &target_rows, &target_cols, &kb_max, &one, temp.getData(), &target_rows, U->getData(), &kb_max, &one, HS->getLocalMatrixElementDown2().getData(), &target_rows);

                                delete imp_kept;

                                delete U;
                                delete UT;
                            }
                        }
                    }
                }
            }
        }
    }
}

void NRG::propagateLocalMatrixElementUp(int n) {
    int Qmax = n + 2;

    //Q loop
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            if(HS != NULL)
            {
                HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q - 1, SZ - 1);
                if(H_other != NULL)
                {
                    //determine dimensions of target matrix
                    int off_target_rows;
                    int off_target_cols;
                    int target_rows;
                    int target_cols;

                    target_rows = HS->getEnergies().getLength();
                    off_target_rows = 0;
                    target_cols = H_other->getEnergies().getLength();
                    off_target_cols = 0;

                    if((target_rows > 0) && (target_cols > 0))
                    {
                        HS->getLocalMatrixElementUp().resize(target_rows, target_cols);
                        HS->getLocalMatrixElementUp().zero();

                        for(int i = 0; i < 4; i++)
                        {
                            int offset_i = 0;
                            int offset_j = 0;
                            for(int g = 0; g < i; g++) {
                                offset_i += HS->getR(g);
                                offset_j += H_other->getR(g);
                            }

                            // variables for matrix multiplication
                            double one = 1.0;
                            double zero = 0.0;
                            char N = 'N';
                            char T = 'T';

                            //3 Matrizes
                            int k_max = HS->getR(i);
                            int kb_max = H_other->getR(i);
                            HilbertSubSpace* H_last = m_hilbertSpaces.getHS(n - 1, Q + dQ[i], SZ + dSz[i]);
                            if((H_last != NULL) && (k_max > 0) && (kb_max > 0))
                            {
                                //copy U matrizes
                                math::Matrix* UT = HS->getHamiltonian().crop(offset_i, off_target_rows, HS->getR(i), target_rows);
                                math::Matrix* U = H_other->getHamiltonian().crop(offset_j, off_target_cols, H_other->getR(i), target_cols);

                                math::Matrix * imp_kept = H_last->getLocalMatrixElementUp().crop(0,0,HS->getR(i), H_other->getR(i));

                                math::Matrix temp;
                                temp.resize(target_rows, kb_max );

                                double factor = signLME[i];

                                lapack::dgemm_(&T, &N, &target_rows, &kb_max, &k_max, &factor, UT->getData(), &k_max, imp_kept->getData() , &k_max, &zero, temp.getData(), &target_rows);
                                lapack::dgemm_(&N, &N, &target_rows, &target_cols, &kb_max, &one, temp.getData(), &target_rows, U->getData(), &kb_max, &one, HS->getLocalMatrixElementUp().getData(), &target_rows);

                                delete imp_kept;

                                delete U;
                                delete UT;
                            }
                        }
                    }
                }
            }
        }
    }
}

void NRG::createPolesF_Up(int n) {
    // within the Andres & Schiller basis it makes only sense to
    // calculate poles for iterations, where states have been already discarded.
    if(n < m_nFirstTruncated)
        return;
    int Qmax = n + 2;
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q-1, SZ - 1);
            if((HS != NULL) && (H_other != NULL))
            {
                math::Matrix& LME = HS->getLocalMatrixElementUp();
                math::Matrix& LME2 = HS->getLocalMatrixElementUp2();
                math::Matrix C(LME);
                int Crows = C.getRows();
                int Ccols = C.getColumns();
#ifdef DEBUG
                if( (Crows < HS->getKeptStates()) || (Ccols < H_other->getKeptStates() ) )
                {
                    std::cout << "ERROR" << std::endl;
                    std::cout << HS->getEnergies().getLength() << std::endl;
                    std::cout << H_other->getEnergies().getLength() << std::endl;
                    throw std::exception();
                }
#endif
                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        C.set(i,j,0.0);

                double zero = 0.0;
                double one = 1.0;
                char N = 'N';
                char T = 'T';

                int m = HS->getEnergies().getLength();
                int n = H_other->getEnergies().getLength();
                int k = H_other->getEnergies().getLength();

                math::Matrix temp;
                temp.resize(m,n);

                lapack::dgemm_(&N, &N, &m, &n, &k, &one, C.getData(), &m, H_other->getDensityMatrixEigenBasis().getData() , &k, &zero, temp.getData(), &m);

                for(int i = 0; i < temp.getRows(); i++)
                {
                    for(int j = 0; j < temp.getColumns(); j++)
                    {
                        double weight = temp.get(i,j) * LME2.get(i,j);
                        double omega = HS->getEnergies().get(i) - H_other->getEnergies().get(j);
                        m_F_Up->addExcitation(omega,weight);
                        //if(calculateExpectationValue)
                        //    expectationValue += weight / (1.0 + exp(omega/m_temperature));
                    }
                }

                math::Matrix temp2;
                temp2.resize(n,m);

                //copy matrix elements
                math::Matrix D(LME2);

                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        D.set(i,j,0.0);

                lapack::dgemm_(&T, &N, &n, &m, &m, &one, D.getData(), &m, HS->getDensityMatrixEigenBasis().getData() , &m, &zero, temp2.getData(), &n);
                for(int i = 0; i < temp2.getRows(); i++)
                {
                    for(int j = 0; j < temp2.getColumns(); j++)
                    {

                        double weight = temp2.get(i,j) * LME.get(j,i);
                        double omega = H_other->getEnergies().get(i) - HS->getEnergies().get(j);
                        m_F_Up->addExcitation(-omega,weight);
                        //if(calculateExpectationValue)
                        //    expectationValue += weight / (1.0 + exp(-omega/m_temperature));
                    }
                }
            }
        }
    }
}

void NRG::createPolesF_Down(int n) {
    // within the Andres & Schiller basis it makes only sense to
    // calculate poles for iterations, where states have been already discarded.
    if(n < m_nFirstTruncated)
        return;
    int Qmax = n + 2;
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q-1, SZ + 1);
            if((HS != NULL) && (H_other != NULL))
            {
                math::Matrix& LME = HS->getLocalMatrixElementDown();
                math::Matrix& LME2 = HS->getLocalMatrixElementDown2();

                math::Matrix C(LME2);
                int Crows = C.getRows();
                int Ccols = C.getColumns();
#ifdef DEBUG
                if( (Crows < HS->getKeptStates()) || (Ccols < H_other->getKeptStates() ) )
                {
                    std::cout << "ERROR" << std::endl;
                    std::cout << HS->getEnergies().getLength() << std::endl;
                    std::cout << H_other->getEnergies().getLength() << std::endl;
                    throw std::exception();
                }
#endif
                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        C.set(i,j,0.0);

                double zero = 0.0;
                double one = 1.0;
                char N = 'N';
                char T = 'T';

                int m = HS->getEnergies().getLength();
                int n = H_other->getEnergies().getLength();
                int k = H_other->getEnergies().getLength();

                math::Matrix temp;
                temp.resize(m,n);

                lapack::dgemm_(&N, &N, &m, &n, &k, &one, C.getData(), &m, H_other->getDensityMatrixEigenBasis().getData() , &k, &zero, temp.getData(), &m);

                for(int i = 0; i < temp.getRows(); i++)
                {
                    for(int j = 0; j < temp.getColumns(); j++)
                    {
                        double weight = temp.get(i,j) * LME.get(i,j);
                        double omega = HS->getEnergies().get(i) - H_other->getEnergies().get(j);
                        m_F_Down->addExcitation(omega,weight);
                        //if(calculateExpectationValue)
                        //    expectationValue += weight / (1.0 + exp(omega/m_temperature));
                    }
                }

                math::Matrix temp2;
                temp2.resize(n,m);

                //copy matrix elements
                math::Matrix D(LME);

                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        D.set(i,j,0.0);

                lapack::dgemm_(&T, &N, &n, &m, &m, &one, D.getData(), &m, HS->getDensityMatrixEigenBasis().getData() , &m, &zero, temp2.getData(), &n);
                for(int i = 0; i < temp2.getRows(); i++)
                {
                    for(int j = 0; j < temp2.getColumns(); j++)
                    {

                        double weight = temp2.get(i,j) * LME2.get(j,i);
                        double omega = H_other->getEnergies().get(i) - HS->getEnergies().get(j);
                        m_F_Down->addExcitation(-omega,weight);
                        //if(calculateExpectationValue)
                        //    expectationValue += weight / (1.0 + exp(-omega/m_temperature));
                    }
                }
            }
        }
    }
}

void NRG::createPolesG_Up(int n) {
    // within the Andres & Schiller basis it makes only sense to
    // calculate poles for iterations, where states have been already discarded.
    if(n < m_nFirstTruncated)
        return;
    int Qmax = n + 2;
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q-1, SZ - 1);
            if((HS != NULL) && (H_other != NULL))
            {
                math::Matrix& LME = HS->getLocalMatrixElementUp();
                math::Matrix& LME2 = HS->getLocalMatrixElementUp();
                math::Matrix C(LME2);
                int Crows = C.getRows();
                int Ccols = C.getColumns();
#ifdef DEBUG
                if( (Crows < HS->getKeptStates()) || (Ccols < H_other->getKeptStates() ) )
                {
                    std::cout << "ERROR" << std::endl;
                    std::cout << HS->getEnergies().getLength() << std::endl;
                    std::cout << H_other->getEnergies().getLength() << std::endl;
                    throw std::exception();
                }
#endif
                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        C.set(i,j,0.0);

                double zero = 0.0;
                double one = 1.0;
                char N = 'N';
                char T = 'T';

                int m = HS->getEnergies().getLength();
                int n = H_other->getEnergies().getLength();
                int k = H_other->getEnergies().getLength();

                math::Matrix temp;
                temp.resize(m,n);

                lapack::dgemm_(&N, &N, &m, &n, &k, &one, C.getData(), &m, H_other->getDensityMatrixEigenBasis().getData() , &k, &zero, temp.getData(), &m);

                for(int i = 0; i < temp.getRows(); i++)
                {
                    for(int j = 0; j < temp.getColumns(); j++)
                    {
                        double weight = temp.get(i,j) * LME.get(i,j);
                        double omega = HS->getEnergies().get(i) - H_other->getEnergies().get(j);
                        m_G_Up->addExcitation(omega,weight);
                        m_n_Up += weight / (1.0 + exp(omega/m_temperature));
                    }
                }

                math::Matrix temp2;
                temp2.resize(n,m);

                //copy matrix elements
                math::Matrix D(LME);

                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        D.set(i,j,0.0);

                lapack::dgemm_(&T, &N, &n, &m, &m, &one, D.getData(), &m, HS->getDensityMatrixEigenBasis().getData() , &m, &zero, temp2.getData(), &n);
                for(int i = 0; i < temp2.getRows(); i++)
                {
                    for(int j = 0; j < temp2.getColumns(); j++)
                    {

                        double weight = temp2.get(i,j) * LME2.get(j,i);
                        double omega = H_other->getEnergies().get(i) - HS->getEnergies().get(j);
                        m_G_Up->addExcitation(-omega,weight);
                        m_n_Up += weight / (1.0 + exp(-omega/m_temperature));
                    }
                }
            }
        }
    }
}

void NRG::createPolesG_Down(int n) {
    // within the Andres & Schiller basis it makes only sense to
    // calculate poles for iterations, where states have been already discarded.
    if(n < m_nFirstTruncated)
        return;
    int Qmax = n + 2;
    for(int Q = Qmax; Q >= -Qmax; Q--)
    {
        int Szmax = (Qmax - abs(Q));
        for(int SZ = -Szmax; SZ <= Szmax; SZ++)
        {
            HilbertSubSpace* HS = m_hilbertSpaces.getHS(n, Q, SZ);
            HilbertSubSpace* H_other = m_hilbertSpaces.getHS(n, Q-1, SZ + 1);
            if((HS != NULL) && (H_other != NULL))
            {
                math::Matrix& LME = HS->getLocalMatrixElementDown();
                math::Matrix C(LME);
                int Crows = C.getRows();
                int Ccols = C.getColumns();
#ifdef DEBUG
                if( (Crows < HS->getKeptStates()) || (Ccols < H_other->getKeptStates() ) )
                {
                    std::cout << "ERROR" << std::endl;
                    std::cout << HS->getEnergies().getLength() << std::endl;
                    std::cout << H_other->getEnergies().getLength() << std::endl;
                    throw std::exception();
                }
#endif
                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        C.set(i,j,0.0);

                double zero = 0.0;
                double one = 1.0;
                char N = 'N';
                char T = 'T';

                int m = HS->getEnergies().getLength();
                int n = H_other->getEnergies().getLength();
                int k = H_other->getEnergies().getLength();

                math::Matrix temp;
                temp.resize(m,n);

                lapack::dgemm_(&N, &N, &m, &n, &k, &one, C.getData(), &m, H_other->getDensityMatrixEigenBasis().getData() , &k, &zero, temp.getData(), &m);

                for(int i = 0; i < temp.getRows(); i++)
                {
                    for(int j = 0; j < temp.getColumns(); j++)
                    {
                        double weight = temp.get(i,j) * LME.get(i,j);
                        double omega = HS->getEnergies().get(i) - H_other->getEnergies().get(j);
                        m_G_Down->addExcitation(omega,weight);
                        m_n_Down += weight / (1.0 + exp(omega/m_temperature));
                    }
                }

                math::Matrix temp2;
                temp2.resize(n,m);

                //copy matrix elements
                math::Matrix D(LME);

                //set kept-kept part to zero
                for(int i = 0; i < HS->getKeptStates(); i++)
                    for(int j = 0; j < H_other->getKeptStates(); j++)
                        D.set(i,j,0.0);

                lapack::dgemm_(&T, &N, &n, &m, &m, &one, D.getData(), &m, HS->getDensityMatrixEigenBasis().getData() , &m, &zero, temp2.getData(), &n);
                for(int i = 0; i < temp2.getRows(); i++)
                {
                    for(int j = 0; j < temp2.getColumns(); j++)
                    {

                        double weight = temp2.get(i,j) * LME.get(j,i);
                        double omega = H_other->getEnergies().get(i) - HS->getEnergies().get(j);
                        m_G_Down->addExcitation(-omega,weight);
                        m_n_Down += weight / (1.0 + exp(-omega/m_temperature));
                    }
                }
            }
        }
    }
}

void NRG::solve_symmetric_SZ(bool silent) {
#ifdef DEBUG
    std::ofstream energiesFile("debug_energies.txt");
#endif
    for(int n = 0; n < m_maxIterations; n++) {
        setupHamiltonian(n);
        truncateStates(n);
        if(n != m_maxIterations-1) {
            propagateChainOperatorElementsUp(n);
            propagateChainOperatorElementsDown(n);
        }
#ifdef DEBUG
        if(n % 2 == 1) {
            energiesFile << n << " ";
            for(int i = 0; i < 10; i++)
                energiesFile << m_energies[i] << " ";
            energiesFile << std::endl;
        }
#endif
    }
    deleteChainOperatorElements(m_maxIterations - 1);
#ifdef DEBUG
    energiesFile.close();
#endif

    // build the density matrix
    builDM();

    // initialize occupations with 0
    m_n_Up = 0.0;
    m_n_Down = 0.0;

    //Second forward run - Propagation of the impurity matrix elements
    for(int n = 0; n < m_maxIterations; n++)
    {
        propagateLocalMatrixElementUp(n);
        propagateLocalMatrixElementUp2(n);
        deleteTransformationMatrices(n);
        createPolesG_Up(n);
        createPolesF_Up(n);
        deleteImpurityMatrixElements(n-1);
        deleteDensityMatrices(n);
    }
    if(!silent)
        std::cout << "n_Up = " << m_n_Up << std::endl;
    m_n_Down = m_n_Up;

    m_F_G_Up = m_G_Up->broaden();
    m_F_F_Up =  m_F_Up->broaden();

    if(!silent) {
        m_G_Up->showStatistics();
        m_F_Up->showStatistics();
    }
    m_F_S_Up = m_U * (m_F_F_Up / m_F_G_Up);

    // correct
    double ImS_corr = 0.0;
    for ( int i = 0; i < m_F_S_Up.getSize(); i++) {
        ImS_corr = std::max(m_F_S_Up.getValueImag(i), ImS_corr);
    }
    for ( int i = 0; i < m_F_S_Up.getSize(); i++) {
        m_F_S_Up.setValue(i, m_F_S_Up.getValue(i) - std::complex<double>(0.0,ImS_corr));
    }
    if(!silent)
        std::cout << "Self-Energy correction: " << ImS_corr << std::endl;
    m_F_S_Down = m_F_S_Up;
    m_F_G_Down = m_F_G_Up;
    m_F_F_Down = m_F_F_Up;
}

void NRG::solve(bool silent) {
    if(m_chainProvider.isSZsymmetric()) {
        solve_symmetric_SZ(silent);
        return;
    }
#ifdef DEBUG
    std::ofstream energiesFile("debug_energies.txt");
#endif
    for(int n = 0; n < m_maxIterations; n++) {
        // std::cout << "iteration: " << n << std::endl;
        setupHamiltonian(n);
        truncateStates(n);
        if(n != m_maxIterations-1) {
            propagateChainOperatorElementsUp(n);
            propagateChainOperatorElementsDown(n);
        }
#ifdef DEBUG
        if(n % 2 == 1) {
            energiesFile << n << " ";
            for(int i = 0; i < 10; i++)
                energiesFile << m_energies[i] << " ";
            energiesFile << std::endl;
        }
#endif
    }
    deleteChainOperatorElements(m_maxIterations - 1);
#ifdef DEBUG
    energiesFile.close();
#endif

    // build the density matrix
    builDM();

    // initialize occupations with 0
    m_n_Up = 0.0;
    m_n_Down = 0.0;

    //Second forward run - Propagation of the impurity matrix elements
    for(int n = 0; n < m_maxIterations; n++)
    {
        propagateLocalMatrixElementUp(n);
        propagateLocalMatrixElementDown(n);
        propagateLocalMatrixElementUp2(n);
        propagateLocalMatrixElementDown2(n);
        deleteTransformationMatrices(n);
        createPolesG_Up(n);
        createPolesG_Down(n);
        createPolesF_Up(n);
        createPolesF_Down(n);
        deleteImpurityMatrixElements(n-1);
        deleteDensityMatrices(n);
    }
    if(!silent){
        std::cout << "n_Up = " << m_n_Up << std::endl;
        std::cout << "n_Down = " << m_n_Down << std::endl;
    }
    m_F_G_Up = m_G_Up->broaden();
    m_F_G_Down = m_G_Down->broaden();
    m_F_F_Up =  m_F_Up->broaden();
    m_F_F_Down =  m_F_Down->broaden();
    m_F_S_Up = m_U * (m_F_F_Up / m_F_G_Up);
    m_F_S_Down = m_U * (m_F_F_Down / m_F_G_Down);

    // correct
    double ImS_corr = 0.0;
    for ( int i = 0; i < m_F_S_Up.getSize(); i++) {
        ImS_corr = std::max(m_F_S_Up.getValueImag(i), ImS_corr);
        ImS_corr = std::max(m_F_S_Down.getValueImag(i), ImS_corr);
    }
    for ( int i = 0; i < m_F_S_Up.getSize(); i++) {
        m_F_S_Up.setValue(i, m_F_S_Up.getValue(i) - std::complex<double>(0.0,ImS_corr));
        m_F_S_Down.setValue(i, m_F_S_Down.getValue(i) - std::complex<double>(0.0,ImS_corr));
    }
    //    SUp.write("sUp.dat");
    //    SDown.write("sDown.dat");

    if(!silent)
        std::cout << "Self-Energy correction: " << ImS_corr << std::endl;

    /* m_G_Up->showInfo();
    math::CFunction GUp = m_G_Up->broaden();
    GUp.write("gUp.dat");
    m_G_Down->showInfo();
    math::CFunction GDown = m_G_Down->broaden();
    GDown.write("gDown.dat");
    m_F_Up->showInfo();
    math::CFunction FUp =  m_F_Up->broaden();
    FUp.write("fUp.dat");
    m_F_Down->showInfo();
    math::CFunction FDown =  m_F_Down->broaden();
    FDown.write("fDown.dat");
    math::CFunction SUp = m_U * (FUp / GUp);
    math::CFunction SDown = m_U * (FDown / GDown);
    SUp.write("sUp.dat");
    SDown.write("sDown.dat");
*/
}

void NRG::showInfo() {
    m_chainProvider.showInfo();

    std::cout << "==========" << std::endl;
    std::cout << "** NRG: **" << std::endl;
    std::cout << "==========" << std::endl;
    std::cout << std::endl;

    std::cout << " - physical parameters:" << std::endl;
    std::cout << " * eps_F = " << m_epsF << " (on-site energy)" << std::endl;
    std::cout << " * U = " << m_U << " (on-site interaction)" << std::endl;
    std::cout << " * T = " << m_temperature << " (temperature)" << std::endl;
    std::cout << " * max. iterations = " << m_maxIterations << " (max. iterations)" << std::endl;
    std::cout << std::endl;
    std::cout << " - numerical parameters:" << std::endl;
    std::cout << " * clusterEnergy = " << m_clusterEnergy << " (cluster energy)" << std::endl;
    std::cout << " * energy cut-off = " << m_energyCutOff << " (energy cut-off)" << std::endl;
    std::cout << " * max. HS dimensions = " << m_maxHSdimension << " (maximum size of Hilbert space)" << std::endl;
    std::cout << " * Enforcing PH symmetry = " << (m_chainProvider.isPHsymmetric() ? "true" : "false") << std::endl;
    std::cout << " * Enforcing SZ symmetry = " << (m_chainProvider.isSZsymmetric() ? "true" : "false") << std::endl;
    std::cout << std::endl;

    m_broadener.showInfo();

}

NRG::~NRG() {
    if(m_G_Up != 0)
        delete m_G_Up;
    if(m_G_Down != 0)
        delete m_G_Down;
    if(m_F_Up != 0)
        delete m_F_Up;
    if(m_F_Down != 0)
        delete m_F_Down;
}

}
