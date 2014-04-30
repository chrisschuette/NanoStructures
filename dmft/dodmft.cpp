#include "dodmft.h"
#include "../config/configuration.h"
#include "dmft.h"

#include <iostream>
#include <memory>

namespace dmft {
void doDMFT(boost::program_options::variables_map vm) {
    std::string configFile;
    config::Configuration& config = config::Configuration::getInstance();
    if(vm.count("config")) {
        configFile = vm["config"].as<std::string>();
        config.readFile(configFile);
    }


    dmft::DMFT DMFT;
    if(vm.count("config"))
        DMFT.configure(config);

    if(vm.count("scaling"))
        DMFT.setScaling(vm["scaling"].as<double>());
    if(vm.count("initialS")) {
        DMFT.setStartWithDelta(false);
        math::CFunction S;
        S.readBinary(vm["initialS"].as<std::string>());
        DMFT.setSelfEnergy(S);
        double scaling = -S.getArgument(0);
        DMFT.setScaling(scaling);
    } else if(vm.count("initialD")) {
        DMFT.setStartWithDelta(true);
        math::CFunction D;
        D.readBinary(vm["initialD"].as<std::string>());
        DMFT.setHybridizationFunction(D);
    }
    if(vm.count("T"))
      DMFT.setT(vm["T"].as<double>());
    if(vm.count("U"))
        DMFT.setU(vm["U"].as<double>());
    if(vm.count("mu"))
        DMFT.setMu(vm["mu"].as<double>());
    if(vm.count("temperature"))
        DMFT.setT(vm["temperature"].as<double>());
    if(vm.count("delta"))
        DMFT.setDelta(vm["delta"].as<double>());
    if(vm.count("tolerance"))
        DMFT.setTolerance(vm["tolerance"].as<double>());
    if(vm.count("maxIterations"))
        DMFT.setMaxIterations(vm["maxIterations"].as<int>());

    std::string outputDir = "";
    if(vm.count("output")) {
        outputDir = vm["output"].as<std::string>();
        if(outputDir[outputDir.size()-1] != '/') {
            outputDir.append("/");
        }
        std::cout << "output directory: " << outputDir << std::endl;
    }

    DMFT.solve(outputDir);

    //output
    DMFT.getGreensFunction().write(outputDir + "G.txt");
    DMFT.getGreensFunction().writeBinary(outputDir + "G.dat");
    DMFT.getSelfEnergy().write(outputDir + "S.txt");
    DMFT.getSelfEnergy().writeBinary(outputDir + "S.dat");
    DMFT.getHybridizationFunction().write(outputDir + "D.txt");
    DMFT.getHybridizationFunction().writeBinary(outputDir + "D.dat");
}
}
