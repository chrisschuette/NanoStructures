#include "hilbertspacetable.h"
#include "hilbertsubspace.h"

namespace nrg {
HilbertSpaceTable::HilbertSpaceTable()
    : m_hilbertSpaces(new HilbertSubSpace* [(MAXITER+1) * (2 * MAXQ + 1) * (2 * MAXSZ + 1)])
{
    for(int i = 0; i < ((MAXITER+1) * (2 * MAXQ + 1) * (2 * MAXSZ + 1)); i++)
        m_hilbertSpaces[i] = NULL;
}

HilbertSpaceTable::~HilbertSpaceTable() {
    for(int i = 0; i < ((MAXITER+1) * (2 * MAXQ + 1) * (2 * MAXSZ + 1)); i++)
        if(m_hilbertSpaces[i] != NULL)
            delete m_hilbertSpaces[i];
    delete [] m_hilbertSpaces;
}

}
