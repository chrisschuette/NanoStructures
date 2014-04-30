#ifndef HYBRIDIZATIONPROVIDER_H
#define HYBRIDIZATIONPROVIDER_H

#include "../chainprovider.h"
#include "../../math/pfunction.h"
#include "../../err/exception.h"

namespace math {
class CFunction;
}

namespace nrg {
namespace chain {
class HybridizationProvider : public ChainProvider
{
public:
    HybridizationProvider();
    void setHybridization(const math::CFunction& hUp, const math::CFunction& hDown);
    void setHybridization(const math::CFunction& hUp);
    math::PFunction getHybridizationUp() { return m_hybridizationUp; }
    math::PFunction getHybridizationDown() { return m_hybridizationDown; }
    virtual void buildChain();
    virtual ~HybridizationProvider();
protected:
    tWilsonChain buildChain(math::PFunction &h);
    mpfr::mpreal sign(mpfr::mpreal a, mpfr::mpreal b);

    math::PFunction m_hybridizationUp;
    math::PFunction m_hybridizationDown;
};

class HybridizationException : error::Exception
{
public:
    HybridizationException(std::string msg) : error::Exception(msg) {}
    virtual ~HybridizationException() throw() {}
};
}
}
#endif // HYBRIDIZATIONPROVIDER_H
