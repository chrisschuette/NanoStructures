#include "flatband.h"

#include <cmath>
#include <iostream>

namespace nrg {
namespace chain {
FlatBand::FlatBand()
    : nrg::ChainProvider()
    , m_spinPolarization(0)
    , m_V(0.1)
{
}

void FlatBand::setSpinPolarization(double polarization) {
    if((polarization != 0.0) && m_symmetry_SZ) {
        std::cerr << "Can't set spin polarization but have SZ symmetry enforced." << std::endl;
        throw std::exception();
    }
    m_spinPolarization = polarization;
}


void FlatBand::buildChain() {
    double V_up = m_V * std::sqrt( 0.5 * ( 1 + m_spinPolarization ) );
    double V_down = m_V * std::sqrt( 0.5 * ( 1 - m_spinPolarization ) );

    for (int n = -1; n < m_length; n++) //do n=-1,m_length-1 = m_length + 1 elements
    {
        if(n==-1)
        {
            m_up.push_back( std::pair< double, double >( V_up, 0.0 ) );
            m_down.push_back( std::pair< double, double >( V_down, 0.0 ) );
        }
        else
        {
            double t = (1.0 - pow(m_lambda,(double) (-n-1) )) /
            ( sqrt( 1.0 - pow(m_lambda,-2*n-1)) * sqrt( 1.0 - pow(m_lambda,-2*n-3) ) );
            m_up.push_back( std::pair< double, double >( t, 0.0 ) );
            m_down.push_back( std::pair< double, double >( t, 0.0 ) );
        }
    }

}

void FlatBand::showInfo() {
    ChainProvider::showInfo();
    std::cout << " - physical parameters:" << std::endl;
    std::cout << " * hybridizationV = " << m_V << " (hybridization strength)" << std::endl;
    std::cout << " * spinPolarization = " << m_spinPolarization << " (spin polarization)" << std::endl;
    std::cout << std::endl;

}



}

}
