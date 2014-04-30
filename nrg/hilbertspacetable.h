#ifndef HILBERTSPACETABLE_H
#define HILBERTSPACETABLE_H

#include <assert.h>
#include <stdlib.h>     /* abs */

#define MAXITER 202
#define MAXQ 20
#define MAXSZ 20

namespace nrg {
class HilbertSubSpace;
class HilbertSpaceTable
{
public:
    HilbertSpaceTable();
    inline HilbertSubSpace* getHS(int iteration, int Q, int Sz) {
      if((abs(Q) <= MAXQ) && (abs(Sz) <= MAXSZ) && (iteration < MAXITER) && (iteration >= -1))
            return m_hilbertSpaces[ iteration + 1 + (Q+MAXQ)*(MAXITER+1)+ (Sz+MAXSZ)*(MAXITER+1)*(2 * MAXQ + 1) ];
        else
            return NULL;
    }
    inline void setHS(int iteration, int Q, int Sz, HilbertSubSpace* HS) {
        assert(abs(Q) <= MAXQ);
        assert(abs(Sz) <= MAXSZ);
        assert(iteration < MAXITER);
        assert(iteration >= -1);
        m_hilbertSpaces[ iteration + 1 + (Q+MAXQ)*(MAXITER+1)+ (Sz+MAXSZ)*(MAXITER+1)*(2 * MAXQ + 1) ] = HS;
    }
    ~HilbertSpaceTable();

protected:
    HilbertSubSpace** m_hilbertSpaces;
};
}

#endif // HILBERTSPACETABLE_H
