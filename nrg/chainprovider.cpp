#include "chainprovider.h"

#include <iostream>

namespace nrg {
ChainProvider::ChainProvider()
    : m_symmetry_PH(true)
    , m_symmetry_SZ(true)
    , m_lambda(2.5)
    , m_length(0)
{
}

void ChainProvider::configure(config::Configuration& configuration) {
    try {
        m_lambda = configuration.getDouble("ChainProvider.lambda");
    } catch( libconfig::SettingNotFoundException ) {}
    try {
        m_length = configuration.getDouble("ChainProvider.chainLength");
    } catch( libconfig::SettingNotFoundException ) {}
    try {
        m_symmetry_PH = configuration.getBool("ChainProvider.PH");
    } catch( libconfig::SettingNotFoundException ) {}
    try {
        m_symmetry_SZ = configuration.getBool("ChainProvider.SZ");
    } catch( libconfig::SettingNotFoundException ) {}
}

void ChainProvider::dump() {
    for(unsigned int i = 0; i < m_up.size(); i++) {
        std::cout << i << " " << m_up.at(i).first << " " << m_up.at(i).second << " " << m_down.at(i).first << " " << m_down.at(i).second << std::endl;
    }
}

void ChainProvider::showInfo() {
    std::cout << "====================" << std::endl;
    std::cout << "** ChainProvider: **" << std::endl;
    std::cout << "====================" << std::endl;
    std::cout << std::endl;

    std::cout << " - numerical parameters:" << std::endl;
    std::cout << " * lambda = " << m_lambda << " (logarithmic discretization)" << std::endl;
    std::cout << " * length = " << m_length << " (wilson chain length)" << std::endl;
    std::cout << std::endl;
}

}
