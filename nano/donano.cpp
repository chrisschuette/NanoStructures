#include "donano.h"

#include "../config/configuration.h"
#include "nanostructure.h"
#include "../mpi/openmpi.h"
#include <memory>
#include <iostream>

namespace nano {
void doNano(boost::program_options::variables_map vm) {

    std::auto_ptr<nano::NanoStructure> structure;
    if(vm.count("structure"))
        structure = std::auto_ptr<NanoStructure>(new nano::NanoStructure(vm["structure"].as<std::string>()));
    else
        structure = std::auto_ptr<NanoStructure>(new nano::NanoStructure);

    std::string configFile;
    config::Configuration& config = config::Configuration::getInstance();
    if(vm.count("config")) {
        configFile = vm["config"].as<std::string>();
        config.readFile(configFile);
	std::cout << "Configuring structure..." << std::endl;
        structure->configure(config);
    }


    if(vm.count("calcCond"))
        structure->setCalculateRhoXX(true);
    if(vm.count("calcHall"))
        structure->setCalculateRhoXY(true);
    if(vm.count("verbose"))
        structure->setVerbose(true);
    if(vm.count("symmetricNano"))
        structure->setSymmetricNano(true);
    // ...

    if(vm.count("resume"))
        structure->readBinary(vm["resume"].as<std::string>());

    if(vm.count("alphaV"))
        structure->setAlphaV(vm["alphaV"].as<double>());

    if(vm.count("maxIterations"))
        structure->setMaxIterations(vm["maxIterations"].as<int>());
    if(vm.count("doECR"))
        structure->setDoECR(true);

    HEAD structure->showInfo();

    if(vm.count("mpi"))
        structure->solve(mpi::OpenMPI::getInstance());
    else
        structure->solve();

}
}
