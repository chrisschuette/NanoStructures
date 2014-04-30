#ifndef CHAINPROVIDER_H
#define CHAINPROVIDER_H

#include "../config/configuration.h"

#include <utility>
#include <vector>

namespace nrg {
namespace chain {
                     // t , // e
typedef std::pair < double, double > tOrbital;
typedef std::vector < tOrbital > tWilsonChain;
}
// epsF   o(0).fs o(1).snd  o(1).fst
// <*>       -       <>       -        <>        -        <>    -    ...
//           V       e0       t0       ...

class ChainProvider
{
public:
    ChainProvider();
    virtual void configure(config::Configuration& configuration);

    virtual void showInfo();

    // getters and setters
    void setLength(int length) { m_length = length; }
    void setLambda(double lambda) { m_lambda = lambda; }
    double getLambda() { return m_lambda; }
    int getLength() { return m_length; }
    inline bool isPHsymmetric() { return m_symmetry_PH; }
    inline bool isSZsymmetric() { return m_symmetry_SZ; }

    inline void setPHsymmetric(bool PH) { m_symmetry_PH = PH; }
    inline void setSZsymmetric(bool SZ) { m_symmetry_SZ = SZ; }

    // retrieve orbital data
    const chain::tOrbital& getOrbitalUp(int n) const { return m_up.at(n); }
    const chain::tOrbital& getOrbitalDown(int n) const { return m_down.at(n); }

    virtual void buildChain() = 0;

    void dump();
protected:
    // particle hole symmetry -- *.second = 0
    bool m_symmetry_PH;
    // sz symmetry -- m_up = m_down
    bool m_symmetry_SZ;

    // energy discretization constant
    double m_lambda;

    // ordered length
    int m_length;

    // chains for spin-up and spin-down band
    chain::tWilsonChain m_up;
    chain::tWilsonChain m_down;
};
}
#endif // CHAINPROVIDER_H
