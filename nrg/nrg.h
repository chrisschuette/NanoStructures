#ifndef NRG_H
#define NRG_H

/**
 * @file
 *
 * @ingroup nrg
 *
 * @brief density matrix implementation of the numerical renormalization group 
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include "../config/configuration.h"
#include "chainprovider.h"
#include "hilbertspacetable.h"
#include "broadener.h"

namespace nrg {
    /*! \brief Density matrix implementation of the numerical renormalization group
     * 
     *  The Numerical Renormalization Group
     *  =================================== 
     *  The Numerical Renormalization Group (NRG) is one of the standard tools
     *  to study correlation effects in quantum impurity models. Here a small
     *  interacting subsystem with a small number of degrees of freedom
     * (the impurity) is coupled to a bath of fermions. No
     * restriction exists as to the structure of the impurity subsystem. 
     * The bath however must consist of non-interacting fermions.
     * 
     * ### Logarithmic discretization ###
     * At the heart of the NRG lies a logarithmic discretization of the continous
     * conduction electron band. The continuum density of states \f$\rho(\epsilon)\f$
     * is approximated by a discrete set of delta poles. By introducing a
     * discretisation parameter \f$\Lambda\f$ Wilson divided the normalised energy
     * range \f$[-1,1]\f$ into \f$2n\f$ intervals where the \f$n\f$th interval
     * (for positive \f$\epsilon\f$) extends from \f$\Lambda^{-(n+1)}\f$ to
     * \f$\Lambda^n\f$. The logarithmic discretisation separates the electron
     * energies into different orders of magnitude where energies close to the
     * Fermi level \f$k_B T \ll D\f$ with \f$D\f$ the bandwidth, which determine
     * the low temperature properties, are well sampled.
     * 
     * ### Mapping onto a chain ###
     * In order to solve the discretized Hamiltonian iteratively one introduces
     * a set of operators \f$f_{n\sigma}\f$ with \f$n>0\f$ in such a way that they exhibit
     * only nearest neighbour coupling. The Hamiltonian assumes the structure of a
     * hopping Hamiltonian on a semi-infinite chain, which is often referred to
     * as the *Wilson chain*.
     * 
     * ### Iterative diagonalization
     * The transformations performed so far have rendered a form of the
     * Hamiltonian which is amendable to an iterative diagonalisation procedure.
     * Now one defines a sequence of Hamiltonians \f$H_N\f$ with \f$N \geq 0\f$.
     *  The full discrete Hamiltonian is recovered in the limit \f$N\to\infty\f$ as 
     * 
     *  \f$ \begin{equation}
     *      H = \lim_{N\to\infty} \frac{1}{2} (1+\Lambda^{-1}) D \Lambda^{-(N-1)/2}H_N
     *      \end{equation}
     * \f$
     * 
     * Two successive Hamiltonians in the series are connected by the recursion relation
     *
     *  \f$
     *      \begin{equation}
     *      H_{N+1} = \Lambda^{1/2} H_N + t_N ( f^\dagger_{N \sigma} f_{N+1\sigma} + f^\dagger _{N+1\sigma} f_{N\sigma} )
     *      \end{equation}
     * \f$
     * 
     * with the initial Hamiltonian in the series containing the impurity itself
     * given by
     * 
     * \f$
     *      \begin{equation}
     *      H_0 = \Lambda^{-\frac{1}{2}} \left[ \tilde{\delta_d} c^\dagger_{d \sigma} c_{d \sigma}  + \tilde{\Gamma}^{1/2} ( f^\dagger_{0 \sigma} c_{d\sigma} + \text{h.c.} ) + \tilde{U} ( c_{d \sigma}^\dagger c_{d \sigma} - 1)^2\right]
     *      \end{equation}
     * \f$
     * 
     * In this form the single impurity Anderson model can be efficiently solved
     * on a computer by taking advantage of the renormalisation group character
     * of the above description. One starts with a diagonalisation of \f$H_0\f$
     * which can be easily accomplished numerically. Assuming that we have
     * diagonalised a Wilson chain of length \f$m\f$ and that the eigenstates are
     * given by \f$|\mathbf r;m\rangle\f$ we construct a product basis for the Wilson
     * chain of length \f$m+1\f$ by
     * 
     * \f$
     *      \begin{equation}
     *      |(\mathbf r,\alpha_{m+1});m+1 \rangle = |\mathbf r; m \rangle  \otimes  |\alpha_{m+1}\rangle
     *      \end{equation}
     * \f$
     * 
     * where \f$|\alpha_{m+1}\rangle\f$ are the eigenstates of the decoupled site
     * \f$|\alpha_{m+1}\rangle=\{ |\rangle,|\uparrow\rangle,|\downarrow\rangle,| \uparrow \downarrow \rangle \}\f$.
     * The matrix elements of the Hamiltonian for the Wilson chain of length \f$m+1\f$
     * for this product basis are given by
     * 
     * \f$
     *   \begin{align}
     *   \langle(\mathbf r',\alpha'_{m+1}); m+1|&H_{m+1} |(\mathbf r,\alpha_{m+1}); m+1\rangle = \Lambda^{1/2} E_{\mathbf r,m} \delta_{\mathbf r \mathbf r'} \delta_{\alpha \alpha'} \nonumber \\
     *   &+\left( \langle \mathbf r'; m |f^\dagger_{m\sigma} |\mathbf r;m \rangle \langle \alpha' | f_{m+1\alpha} |\alpha \rangle + \langle \mathbf r'; m | f_{m\sigma} |\mathbf r;m \rangle \langle\alpha'| f^\dagger_{m+1\alpha} |\alpha\rangle\right)
     *   \end{align}
     * \f$
     * 
     * The eigenvalue problem for the chain of length \f$m+1\f$ can therefore be
     * solved numerically using only a knowledge of the eigenspectrum of the
     * chain of length \f$m\f$ and the matrix elements of the operators \f$f_{m\sigma}^\dagger\f$.
     * Diagonalising the Hamiltonian \f$H_{m+1}\f$, set up in the above product basis,
     * can be described by a unitary transformation
     * 
     * \f$
     *   \begin{equation}
     *   |\mathbf r'; m+1 \rangle = \sum_{\alpha_{m+1}, \mathbf r} U^{\alpha_{m+1}}_{\mathbf r', \mathbf r} |\mathbf r;m\rangle \otimes |\alpha_{m+1}\rangle
     *   \end{equation}
     * \f$
     * 
     * where \f$|\mathbf r'; m+1 \rangle\f$ denotes the new eigenbasis of the
     * Hamiltonian \f$H_{m+1}\f$.
     * 
     * \image html truncation.png "Illustration of the truncation procedure. The
     * iterative diagonalisation splits each energy level into 4 levels upon
     * the addition of another chain element. In this schematic picture however 
     *each energy level is split into only two levels in order not to make the 
     * illustration to cluttered. Due to the exponential decrease in the 
     * couplings it is save to truncate the high energy states without altering 
     * the spectrum of the low energy states. The truncated states are marked 
     * red."
     * 
     * References
     * ----------
     * * <a href="http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.80.395">NRG review article by R. Bulla et al.</a>
     * * <a href="http://journals.aps.org/prb/abstract/10.1103/PhysRevB.21.1003">Renormalization-group approach to the Anderson model  by  H. R. Krishna-murthy, J. W. Wilkins, and K. G. Wilson</a> 
     */
    class NRG {

        enum STATES {
            KEPT, DISCARDED, BOTH
        };

    public:
        /**
         * \brief constructs an NRG algorithm instance
         * The Wilson chain is provided by an implementation of the ChainProvider
         * interface. The Broadener is used after diagonalization of the chain
         * and calculation of the density matrix to broaden the representation
         * of various correlation functions as a discrete set of delta peaks
         * into continous functions.
         * 
         * \param[in] chainProvider implementation of the ChainProvider interface;
         *                          supplies the Wilson chain.
         * \param[in] broadener     implementation of the Broadener interface;
         *                          broadens the discrete set of delta peaks
         *                          of the correlation functions into continous
         *                          functions.
         */
        NRG(ChainProvider& chainProvider, Broadener& broadener);
  
        /**
         * \brief destructs the NRG algorithm object and frees all previously
         *        allocated memory.
         */
        ~NRG();

        
        /**
         * \brief configures the NRG instance with information in the configuration object.
         * @param[in] configuration a configuration object which can be queried
         *                          for various numerical and physical parameters
         *                          concerning the impurity problem.
         */
        void configure(config::Configuration& configuration);

        // setters and getters -- physical parameters

        /**
         * \brief retrieves the strength of the on-site impurity interaction U
         * @return on-site impurity interaction U
         */
        double getU() const {
            return m_U;
        }

        /**
         * \brief sets the strength of the on-site impurity interaction U
         * @param U on-site impurity interaction U
         */
        void setU(double U) {
            m_U = U;
        }

        /**
         * \brief retrieves the on-site impurity energy \f$\epsilon_f\f$
         * @return epsF on-site energy \f$\epsilon_f\f$
         */
        double getEpsF() const {
            return m_epsF;
        }

        /**
         * \brief sets the on-site impurity energy \f$\epsilon_f\f$
         * @param[in] epsF on-site energy \f$\epsilon_f\f$
         */
        void setEpsF(double epsF) {
            m_epsF = epsF;
        }

        /**
         * \brief retrieves the system temperature T
         * @return system temperature T
         */
        double getTemperature() const {
            return m_temperature;
        }

        /**
         * \brief sets the system temperature T
         * @param[in] T system temperature T
         */
        void setTemperature(double T) {
            m_temperature = T;
        }

        //setters and getters -- numerical parameters
        
        /**
         * \brief returns the cluster energy.
         * 
         * Symmetries create degeneracies in the spectrum of eigenstates. Due to numerical noise
         * these states can energetically drift apart by tiny amounts. It is generally
         * not a good idea to cut-off within a multiplet (or "cluster"). It is ensured that no
         * cut is made in between two energy levels less than a certain energy difference
         * apart. This function returns this energy difference.
         * 
         * @return cluster energy
         */
        double getClusterEnergy() const {
            return m_clusterEnergy;
        }
        /**
         * \brief sets the cluster energy to the indicated value. see getClusterEnergy().
         * 
         * @param clusterEnergy
         */
        void setClusterEnergy(double clusterEnergy) {
            m_clusterEnergy = clusterEnergy;
        }

        /**
         * \brief returns the high energy cut-off for state truncation
         * 
         * @return high-energy cut-off
         */
        double getEnergyCutOff() const {
            return m_energyCutOff;
        }

        /**
         * \brief sets the high energy cut-off for state truncation
         * 
         * @param[in] clusterEnergy new high-energy cut-off
         */
        void setEnergyCutOff(double energyCutOff) {
            m_energyCutOff = energyCutOff;
        }

        /**
         * \brief  returns the maximum allowed size of the HilberSpace. See truncateStates(int n).
         * 
         * @return max dimension of HilbertSpace
         */
        int getMaxHilbertSpaceDimension() {
            return m_maxHSdimension;
        }
        
        /**
         * \brief  sets the maximum allowed size of the HilberSpace. See truncateStates(int n).
         * 
         * It is advisable to set this number to a large value (so that it does not
         * get in the way) and let the energy cut-off( see setEnergyCutOff(double energyCutOff) )
         * do the truncation.
         * 
         * @param[in] maxHSdimension max. dimension of the HilbertSpace
         */
        void setMaxHilbertSpaceDimension(int maxHSdimension) {
            m_maxHSdimension = maxHSdimension;
        }

        /**
         * returns the number of NRG iterations to be performed.
         * 
         * @return NRG iterations
         */
        int getMaxIterations() {
            return m_maxIterations;
        }

        /**
         * sets the number of NRG iterations to be performed.
         * 
         * @param[in] maxIterations # of NRG iterations
         */
        void setMaxIterations(int maxIterations) {
            m_maxIterations = maxIterations;
        }

        /**
         * \brief returns the exp. value for spin-\f$\uparrow\f$ electrons
         *        on the impurity.
         * The actual calculation is performed in createPolesG_Up(int n).
         * @return 
         */
        inline double getOccupationUp() {
            return m_n_Up;
        }

        /**
         * \brief returns the exp. value for spin-\f$\downarrow\f$ electrons
         *        on the impurity.
         * The actual calculation is performed in createPolesG_Down(int n).
         * @return 
         */
        inline double getOccupationDown() {
            return m_n_Down;
        }

        /**
         * \brief returns the exp. value for electrons
         *        on the impurity.
         * The actual calculation is performed in createPolesG_Up(int n)
         * and createPolesG_Down(int n).
         * @return 
         */
        inline double getOccupation() {
            return m_n_Up + m_n_Down;
        }

        /**
         * \brief returns the exp. value for the magnetization
         * 
         * The actual calculation is performed in createPolesG_Up(int n)
         * and createPolesG_Down(int n). The magnetization is calculated as
         * \f$ \avg{m} = \avg{n_\uparrow} - \avg{\downarrow}\f$
         * @return 
         */
        inline double getMagnetization() {
            return m_n_Up - m_n_Down;
        }

        // nrg
        /** 
         * \brief prepares the NRG instance for the iterative diagonalization.
         * 
         * The function performs the following steps:
         * #### Initialization ####
         * * If the max. iterations have not been set, chooses them automatically
         *   based on the system temperature T.
         * * Checks whether enough NRG iteration will be performed for the given
         *   temperature T.
         * * Calls the ChainProvider to prepare a chain of the appropriate length,
         * * Sets the temperature T in the broadeners and initializes them.
         * * Calls setupInitialState() to initialize the impurity Hamiltonian \f$H_0\f$.
         */
        void init();

        /**
         * \brief initializes the impurity Hamiltonian \f$H_0\f$.
         * 
         * This function creates the first 4 HilbertSubSpace instances for the chain
         * consisting only of the impurity
         * ( \f$ \{ |\rangle,|\uparrow\rangle,|\downarrow\rangle,| \uparrow \downarrow \rangle \}\f$ ).
         * It determines the eigenenergies and both the local impurity and the
         * chain matrix elements.
         */
        void setupInitialState();
        
        /**
         * \brief iterates over all possible quantum number Q and Sz for iteration N
         *        and constructs the Hamiltonian for the corresponding HilbertSubSpace.
         * 
         * For a given HilbertSubSpace (HSS) in iteration N the function looks up the appropriate HSSs (up to 4)
         * of iteration N-1 and uses the chain matrix elements to construct the Hamiltonian
         * for the HSS in iteration N. Then the Hamiltonian is diagonalized using
         * a call to Lapack's <a href="http://www.netlib.org/lapack/explore-3.1.1-html/dsyev.f.html">dsyev</a>
         * Lapack replaces the Hamiltonian m_H in the HSS by
         * the unitary transformation that diagonalized it and places and sorts
         * the eigenvalues (energies) into the HSS object as well. The function
         * discards the chain matrix elements of iteration N-1, which are no longer needed.
         * 
         * @param[in] iteration NRG iteration
         */
        void setupHamiltonian(int iteration); // iteration >= 0
        
        /**
         * \brief marks all eigenstates above a cut-off energy \f$E_c\f$ as discarded.
         * 
         * The Hilbertspace grows exponentially with the length of the Wilson
         * chain. Therefore the eigenstates have to be truncated above an 
         * energy cutoff (see void setEnergyCutOff(double energyCutOff) ). Only 
         * the "kept" states in iteration N-1 with an energy smaller than the cut-off are used to
         * construct the Hilbertspace for iteration N. A separation of energy scales
         * is achieved by the logarithmic discretization of the energy band and
         * justifies the truncation: States with a high energy in iteration N-1
         * cannot affect the structure of the low-energy states in iteration N.
         * 
         * Special attention has to be paid to the role of symmetries. Symmetries
         * create degeneracies in the spectrum of eigenstates. Due to numerical noise
         * these states can energetically drift apart by tiny amounts. It is generally
         * not a good idea to cut-off within a multiplet. It is ensured that no
         * cut is made in between two energy levels less than a certain energy difference
         * apart (see void setClusterEnergy(double clusterEnergy) ).
         * 
         * @param[in] iteration NRG iteration
         */
        void truncateStates(int iteration);

        
        /**
         * \brief calculates the local impurity matrix elements or spin-\f$\uparrow\f$
         *        for a given iteration.
         * 
         * Lehmann resolving the impurity spectral function one obtains
         * 
         * \f$
         * \begin{equation}
         *  A_\sigma(\omega) = \sum_{a,b} \bra{b} c_{d\sigma}\ket{a} \frac{\exp{\left[ -\beta E_a \right]}}{Z} \bra{a}c_{d\sigma}^\dagger\ket{b} \delta(\omega + E_a - E_b)
         * \end{equation}
         * \f$
         * 
         * where $Z$ is the total partition function
         * $Z = \sum_a \exp{\left[ -\beta E_a \right]}$, $\ket{a}$ and $\ket{b}$
         * are a complete set of states and $E_a$ is the eigenenergy of state
         * $\ket{a}$. We see that the Lehmann representation gathers the
         * necessary information to construct the spectrum from knowledge of
         * certain matrix elements encoding hopping processes between the
         * impurity and the conduction electron band. It is precisely these matrix
         * elements \f$ \langle \mathbf r'; m | c^\dagger_{\uparrow} | \mathbf r;m \rangle \delta_{\alpha' \alpha} \f$
         * which are calculated in this function.
         * 
         * @param[in] iteration NRG iteration
         */
        void propagateLocalMatrixElementUp(int iteration);
        
        /**
         * \brief see propagateLocalMatrixElementUp(int iteration). This is the spin-\f$\downarrow\f$ version.
         * @param[in] iteration NRG iteration
         */
        void propagateLocalMatrixElementDown(int iteration);
        
        /**
         * \brief calculates other local impurity matrix elements or spin-\f$\uparrow\f$
         *        for a given iteration.
         * 
         * Bulla et al. first showed that it is possible to write the self-energy
         * as the ratio of two correlation functions, both of which can be
         * calculated directly within the NRG. An equation of motion technique
         * is used to show that the self-energy is given by,
         * 
         * \f$
         *   \begin{equation}
         *   \Sigma_\sigma(\omega) = U \frac{F_{\sigma}(\omega)}{G_{\sigma}(\omega)}
         *   \end{equation}
         * \f$
         * 
         * where $G_{\sigma}(\omega)$ is the impurity Green's function defined as
         * 
         * \f$
         *   \begin{equation}
         *   G_{\sigma}(\omega) = -i \int_{-\infty}^\infty dt~e^{ i \omega t } \Theta(t) \avg{ \left\{ f_{\sigma}(t), f_{\sigma}^\dagger \right\} }
         *   \end{equation}
         * \f$
         * 
         * and $F_{\sigma}(\omega)$ is an auxiliary correlation function given by
         * 
         * \f$
         *   \begin{equation}
         *   F_{\sigma}(\omega) = -i \int_{-\infty}^\infty dt e^{ i \omega t } \Theta(t) \avg{ \left\{ \left( f_{\bar{\sigma}}^\dagger f_{\bar{\sigma}} f_{\sigma}\right)(t), f_{\sigma}^\dagger \right\} }
         *   \end{equation}
         * \f$
         * 
         * The imaginary parts of $F_{\sigma}(\omega)$ and $G_{\sigma}(\omega)$
         * are calculated from NRG data using the Lehman sum within the full
         * density matrix approach, with the poles of the spectrum broadened as
         * above. The real parts are then obtained by Kramers-Kronig transform.
         * Discretization artifacts cancel to some extent by dividing the two
         * quantities. This produces a rather smooth self-energy, which in term
         * can be used to calculate an improved spectrum for the impurity.
         * Z-averaging can also be used to
         * further increase accuracy and resolution.
         * 
         * It is precisely these matrix
         * elements \f$ \langle \mathbf r'; m | c^\dagger_{\downarrow} c^\phantom{\dagger}_{\downarrow} c^\dagger_{\uparrow} | \mathbf r;m \rangle \delta_{\alpha' \alpha} \f$
         * which are calculated in this function.
         * 
         * @param[in] iteration NRG iteration
         */
        void propagateLocalMatrixElementUp2(int iteration);
        
        /**
         * \brief see propagateLocalMatrixElementUp2(int iteration). This is the spin-\f$\downarrow\f$ version.
         * @param[in] iteration NRG iteration
         */
        void propagateLocalMatrixElementDown2(int iteration);
        
        /**
         * \brief calculates the chain matrix elements for spin-\f$\uparrow\f$
         *        for a given iteration.
         * 
         * Assuming that we have diagonalised a Wilson chain of length \f$m\f$ and
         * that the eigenstates are given by \f$|\mathbf r;m\rangle\f$ we construct a
         * product basis for the Wilson chain of length \f$m+1\f$ by
         * 
         * \f$
            \begin{equation}
            |(\mathbf r,\alpha_{m+1});m+1\rangle = |\mathbf r; m\rangle  \otimes  |\alpha_{m+1}\rangle
            \end{equation}
         * \f$
         * 
         * where \f$|\alpha_{m+1}\rangle\f$ are the eigenstates of the decoupled site
         * \f$|\alpha_{m+1}\rangle=\{ |\rangle,|\uparrow\rangle,|\downarrow\rangle,| \uparrow \downarrow \rangle \}\f$.
         * The matrix elements of the Hamiltonian for the Wilson chain of
         * length \f$m+1\f$ for this product basis are given by
         * 
         * \f$
            \begin{align}
            \langle(\mathbf r',\alpha'_{m+1}); m+1| &H_{m+1} |(\mathbf r,\alpha_{m+1}); m+1\rangle = \Lambda^{1/2} E_{\mathbf r,m} \delta_{\mathbf r \mathbf r'} \delta_{\alpha \alpha'} \nonumber \\
            &+\left( \langle \mathbf r'; m | f^\dagger_{m\sigma} | \mathbf r;m \rangle \langle \alpha'| f_{m+1\alpha} |\alpha \rangle + \langle \mathbf r'; m | f_{m\sigma} |\mathbf r;m\rangle \langle\alpha'| f^\dagger_{m+1\alpha} |\alpha\rangle\right).
            \end{align}
         * \f$
         * 
         * It is precisely the matrix elements \f$ \langle \mathbf r'; m | f^\dagger_{m\uparrow} | \mathbf r;m \rangle \f$,
         * which are calculated in this function.
         *
         * @param[in] iteration NRG iteration
         */
        void propagateChainOperatorElementsUp(int iteration);
        
        /**
         * \brief see propagateChainOperatorElementsUp(int iteration). This is the spin-\f$\downarrow\f$ version.
         * 
         * @param[in] iteration NRG iteration
         */
        void propagateChainOperatorElementsDown(int iteration);

        /**
         * \brief calculates delta peaks for the spin-\f$\uparrow\f$ impurity Green's function.
         * 
         * This function uses the reduced density matrix and the impurity operator matrix
         * elements to calculate delta peaks of the Lehmann representation
         * of the spin-\f$\uparrow\f$ impurity Green's function for the given iteration.
         * 
         * @param[in] iteration NRG iteration
         */
        void createPolesG_Up(int iteration);
        
        /**
         * \brief see createPolesG_Up(int iteration). This is the spin-\f$\downarrow\f$ version.
         * 
         * @param[in] iteration NRG iteration
         */
        void createPolesG_Down(int iteration);
        
        /**
         * \brief calculates delta peaks for the spin-\f$\uparrow\f$ impurity 
         * correlator used in the self-energy trick.
         * 
         * This function uses the reduced density matrix and the impurity operator matrix
         * elements to calculate delta peaks of the Lehmann representation
         * of the spin-\f$\uparrow\f$ correlator for the given iteration used
         * in the self-energy trick.
         *
         * @param[in] iteration NRG iteration
         */
        void createPolesF_Up(int iteration);
        
        /**
         * \brief see createPolesF_Down(int iteration). This is the spin-\f$\downarrow\f$ version.
         * 
         * @param[in] iteration NRG iteration
         */
        void createPolesF_Down(int iteration);

        /**
         * \brief solves the impurity problem.
         * 
         * First, the function queries the ChainProvider whether the provided chain
         * is symmetric with respect to the Sz quantum number. If so, it suffices
         * to calculate all matrix elements and correlator for spin-\f$\uparrow\f$
         * only.
         * 
         * The algorithm moves along the Wilson chain iteratively in a forward-backward-
         * forward pattern.
         * 
         * #### Step 1 - Forward run: Iterative diagonalization ####
         * In the first forward run, the Hamiltonian is iteratively constructed
         * and diagonalized (using setupHamiltonian(int n) ). Chain matrix 
         * elements are calculated along the way as needed. Once the highest
         * eigenenergy of the chain exceeds the predefined threshold value, the truncation
         * sets in (see truncateStates(int iteration) )  and only the kept
         * states are used to construct the HilbertSpace of the next iteration.
         * 
         * #### Step 2 - Backward run: Reduced density matrix ####
         * In a second step, buildDM() is called, which performs a backward run
         * starting from the last iteration and calculates the reduced density 
         * matrix for the different iterations.
         * 
         * #### Step 3 - Forward run: Calculation of Correlators ####
         * In a third step, a forward run is performed which uses the reduced density
         * matrix to calculate impurity matrix elements. The matrices are propagated from
         * iteration to iteration and used to calculate an approximate Lehmann representation
         * of the impurity Green's function and a another correlation function
         * which is needed for the self-energy trick.
         * 
         * Finally, the discrete set of delta peaks for each correlator is broadened
         * to get continous functions (see Broadener).
         * 
         * \param[in] silent if true, silences the impurity solver (useful for DMFT)
         */
        void solve(bool silent = false);
        
        /**
         * \brief Version of solve(bool silent) for a Wilson chain with Sz symmetry.
         * 
         * @param silent if true, silences the impurity solver (useful for DMFT)
         */
        void solve_symmetric_SZ(bool silent = false);

        /**
         * \brief performs a backward run to construct the reduced density matrices
         * 
         */
        void builDM();

        /**
         * \brief lists the physical and numerical parameters for the NRG instance.
         *        Also dumps info about the Broadener object.
         */
        void showInfo();

        // clean-up code
        /**
         * \brief deletes all chain operator matrix elements for the indicated iteration.
         * @param[in] iteration NRG iteration
         */
        void deleteChainOperatorElements(int iteration);
        
        /**
         * \brief deletes all transformation matrices for the indicated iteration.
         * @param[in] iteration NRG iteration
         */        
        void deleteTransformationMatrices(int iteration);
        
         /**
         * \brief deletes all chain operator matrix elements for the indicated iteration.
         * @param[in] iteration NRG iteration
         */
        void deleteDensityMatrices(int iteration);
                
        /**
         * \brief deletes all impurity operator matrix elements for the indicated iteration.
         * @param[in] iteration NRG iteration
         */
        void deleteImpurityMatrixElements(int iteration);

        /**
         * \brief returns the impurity selfenergy.
         * @param[out] SUp reference to spin-\f$\uparrow\f$ self energy
         * @param[out] SDown reference to spin-\f$\downarrow\f$ self energy
         */
        void getSelfEnergy(math::CFunction& SUp, math::CFunction& SDown) {
            SUp = m_F_S_Up;
            SDown = m_F_S_Down;
        }

        /**
         * \brief returns the impurity selfenergy.
         * @param[out] S reference to spin-\f$\uparrow\f$ self energy
         */
        void getSelfEnergy(math::CFunction& S) {
            S = m_F_S_Up;
        }

        /**
         * \brief returns the impurity Green's function.
         * @param[out] GUp reference to spin-\f$\uparrow\f$ Green's function
         * @param[out] GDown reference to spin-\f$\downarrow\f$ Green's function
         */
        void getGreensFunction(math::CFunction& GUp, math::CFunction& GDown) {
            GUp = m_F_G_Up;
            GDown = m_F_G_Down;
        }

        /**
         * \brief returns the Green's function.
         * @param[out] G reference to spin-\f$\uparrow\f$ Green's function
         */
        void getGreensFunction(math::CFunction& G) {
            G = m_F_G_Up;
        }

        /**
         * \brief returns the correlator for the self-energy trick
         * @param FUp reference to spin-\f$\uparrow\f$ correlator
         * @param FDown reference to spin-\f$\downarrow\f$ correlator
         */
        void getFFunction(math::CFunction& FUp, math::CFunction& FDown) {
            FUp = m_F_F_Up;
            FDown = m_F_F_Down;
        }

        /**
         * \brief returns the correlator for the self-energy trick
         * @param F reference to spin-\f$\uparrow\f$ correlator
         */
        void getFFunction(math::CFunction& F) {
            F = m_F_F_Up;
        }

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
        std::mathbfor<double> m_energies;

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
