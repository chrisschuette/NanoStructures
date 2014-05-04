#ifndef HILBERTSPACETABLE_H
#define HILBERTSPACETABLE_H

/**
 * @file
 *
 * @ingroup nrg
 *
 * @brief Hilbert subspace table for charge quantum number Q and spin S
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include <assert.h>
#include <stdlib.h>     /* abs */

#define MAXITER 202
#define MAXQ 20
#define MAXSZ 20

namespace nrg {
    class HilbertSubSpace;

    /*! \brief HilbertSpaceTable: Manages the table of HilbertSubSpace instances
     */
    class HilbertSpaceTable {
    public:
        /**
         * \brief constructs a table of HilbertSubSpace instances.
         * The size of the table is determined by MAXITER, MAXQ, MAXSZ
         */
        HilbertSpaceTable();

        /**
         * \brief returns a pointer to the HilbertSubSpace for a given iteration, charge Q and spin Sz.
         * If the HilbertSubSpace does not exist NULL is returned.
         * \param[in] iteration iteration in the chain diagonalization
         * \param[in] Q charge quantum number
         * \param[in] Sz spin quantum number
         * \return pointer to HilbertSubSpace otherwise NULL
         */
        inline HilbertSubSpace* getHS(int iteration, int Q, int Sz) {
            if ((abs(Q) <= MAXQ) && (abs(Sz) <= MAXSZ) && (iteration < MAXITER) && (iteration >= -1))
                return m_hilbertSpaces[ iteration + 1 + (Q + MAXQ)*(MAXITER + 1)+ (Sz + MAXSZ)*(MAXITER + 1)*(2 * MAXQ + 1) ];
            else
                return NULL;
        }
        
        /**
         * \brief sets the HilbertSubSpace for a given iteration, charge Q and spin Sz.
         * \param[in] iteration iteration in the chain diagonalization
         * \param[in] Q charge quantum number
         * \param[in] Sz spin quantum number
         * \param[in] HS pointer to the HilbertSubSpace
         */
        inline void setHS(int iteration, int Q, int Sz, HilbertSubSpace* HS) {
            assert(abs(Q) <= MAXQ);
            assert(abs(Sz) <= MAXSZ);
            assert(iteration < MAXITER);
            assert(iteration >= -1);
            m_hilbertSpaces[ iteration + 1 + (Q + MAXQ)*(MAXITER + 1)+ (Sz + MAXSZ)*(MAXITER + 1)*(2 * MAXQ + 1) ] = HS;
        }
        
        /**
         * \brief destructs the HilbertSubSpace table and deletes all HilberSubSpace instances.
         */
        ~HilbertSpaceTable();

    protected:
        HilbertSubSpace** m_hilbertSpaces;
    };
}

#endif // HILBERTSPACETABLE_H
