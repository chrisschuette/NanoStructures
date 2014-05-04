#ifndef HILBERTSUBSPACE_H
#define HILBERTSUBSPACE_H

/**
 * @file
 *
 * @ingroup nrg
 *
 * @brief Hilbert subspace for charge quantum number Q and spin S
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include "../math/matrix.h"
#include "../math/vector.h"

namespace nrg {

/*! \brief HilbertSubSpace: Hilbert subspace for charge quantum number Q and spin Sz.
     * This class contains the eigenenergies, matrix elements and reduced density matrix
     * for a Hilbert subspace with charge quantum number Q and spin S.
     */
    class HilbertSubSpace {
    public:
        /**
         * \brief constructs a subspace with quantum numbers Q and Sz
         * \param[in] Q charge quantum number
         * \param[in] Sz spin quantum number
         */
        HilbertSubSpace(int Q, int Sz);

        // getters and setters
        
        /**
         * \brief returns the charge quantum number Q of the subspace
         * \return charge quantum number
         */
        inline int getQ() {
            return m_Q;
        }

        /**
         * \brief returns the spin quantum number Sz of the subspace
         * \return spin quantum number
         */
        inline int getSz() {
            return m_Sz;
        }

        /**
         * \brief returns the numer of kept states for the subspace
         * \return # of kept states
         */
        inline int getKeptStates() {
            return m_nKeptStates;
        }

        /**
         *  \brief sets the number of kept states for the subspace
         * \param[in] keptStates # of kept states
         */
        void setKeptStates(int keptStates) {
            m_nKeptStates = keptSta tes;
        }

        /**
         *  \brief returns reference to the Hamiltonian/unitary transformation
         * After diagonalization Lapack turns this Hamiltonian matrix H into 
         * the diagonalizing unitary transformation U. So beware!
         * \return reference to matrix H/U
         */
        inline math::Matrix& getHamiltonian() {
            return m_H;
        }

        /**
         * \brief set the Hamiltonian H
         * \param[in] H Hamiltonian matrix
         */
        void setHamiltonian(math::Matrix& H) {
            m_H = H;
        }

        /**
         * \brief returns a reference to the vector of eigenenergies
         * \return eigenenergies
         */
        inline math::Vector& getEnergies() {
            return m_E;
        }

        /**
         * \brief sets the # of kept states for ith hilbert space of the prev. iteration
         * \param[in] i 0,1,2,3
         * \param[in] r # of kept states
         */
        void setR(int i, int r) {
            m_r[i] = r;
        }

        /**
         * \brief retrieves the # of kept states for ith hilbert space of the prev. iteration
         * \param[in] i 0,1,2,3
         * \return # of states
         */
        int getR(int i) {
            return m_r[i];
        }

        /**
         * \brief retrieves a reference to chain matrix elements for spin up electrons
         * \return chain matrix elements for spin up
         */
        inline math::Matrix& getChainElementsUp() {
            return m_chainOperatorElementsUp;
        }

        /**
         * \brief retrieves a reference to chain matrix elements for spin down electrons
         * \return chain matrix elements for spin down
         */
        inline math::Matrix& getChainElementsDown() {
            return m_chainOperatorElementsDown;
        }

        inline math::Matrix& getDensityMatrixEigenBasis() {
            return m_densityMatrixEigenBasis;
        }

        inline math::Matrix& getDensityMatrix() {
            return m_densityMatrix;
        }

        /**
         * \brief retrieves a reference to impurity matrix elements for spin down electrons
         * Take a look at the "self-energy trick" by R. Bulla
         * \return impurity matrix elements for spin down
         */
        inline math::Matrix& getLocalMatrixElementUp() {
            return m_fUp;
        }

        /**
         * \brief retrieves a reference to impurity matrix elements for spin down electrons
         * Take a look at the "self-energy trick" by R. Bulla
         * \return impurity matrix elements for spin down
         */
        inline math::Matrix& getLocalMatrixElementDown() {
            return m_fDown;
        }

        /**
         * \brief retrieves a reference to other impurity matrix elements for spin down electrons
         * Take a look at the "self-energy trick" by R. Bulla
         * \return other impurity matrix elements for spin down
         */
        inline math::Matrix& getLocalMatrixElementUp2() {
            return m_fffUp;
        }

        /**
         * \brief retrieves a reference to other impurity matrix elements for spin down electrons
         * Take a look at the "self-energy trick" by R. Bulla
         * \return other impurity matrix elements for spin down
         */
        inline math::Matrix& getLocalMatrixElementDown2() {
            return m_fffDown;
        }
    protected:
        int m_Q;
        int m_Sz;
        int m_r[4];
        int m_nKeptStates;

        math::Matrix m_H;
        math::Vector m_E;
        math::Matrix m_chainOperatorElementsUp;
        math::Matrix m_chainOperatorElementsDown;
        math::Matrix m_densityMatrix;
        math::Matrix m_densityMatrixEigenBasis;
        math::Matrix m_fUp;
        math::Matrix m_fDown;
        math::Matrix m_fffUp;
        math::Matrix m_fffDown;
    };
}

#endif // HILBERTSUBSPACE_H
