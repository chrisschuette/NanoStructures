#ifndef FLATBAND_H
#define FLATBAND_H

#include "../chainprovider.h"

namespace nrg {
namespace chain {
class FlatBand : public ChainProvider
{
public:
    FlatBand();

    virtual void showInfo();


    //getters and setters
    void setSpinPolarization(double polarization);
    double getSpinPolarization() const { return m_spinPolarization; }

    void setHybridization(double V) { m_V = V; }
    double getHybridization() { return m_V; }

    virtual void buildChain();

protected:
    double m_spinPolarization; // -1 ... 0 ... 1
    double m_V;
};
}
}
#endif // FLATBAND_H
