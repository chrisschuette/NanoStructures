#include "hilbertsubspace.h"

namespace nrg {
HilbertSubSpace::HilbertSubSpace(int Q, int Sz)
    : m_Q(Q)
    , m_Sz(Sz)
    , m_nKeptStates(-1) // means all are kept
{
    m_r[0] = 0;
    m_r[1] = 0;
    m_r[2] = 0;
    m_r[3] = 0;
}
}
