#include "donrg.h"

#include "../utils/utils.h"

#include "nrg.h"
#include "broadening/fftbroadener.h"
#include "chain/flatband.h"
#include "chain/hybridizationprovider.h"

#include "../config/configuration.h"

#include <memory>

namespace nrg {
void doNRG(boost::program_options::variables_map vm) {
    std::string configFile;
    config::Configuration& config = config::Configuration::getInstance();
    if(vm.count("config")) {
        std::string configFile = vm["config"].as<std::string>();
        config.readFile(configFile);
    }
    //RAII
    std::auto_ptr<nrg::NRG> NRG;
    std::auto_ptr<nrg::Broadener> broadener;
    std::auto_ptr<nrg::ChainProvider> chain;

    // create & configure chain provider
    if(!vm.count("chain") || (utils::toUpperCase(vm["chain"].as<std::string>()) == "FLAT")) {
        nrg::chain::FlatBand* flatBand = new nrg::chain::FlatBand;
        if(vm.count("config"))
            flatBand->configure(config);
        if(vm.count("PHsymmetry"))
            flatBand->setLambda(vm["PHsymmetry"].as<bool>());
        if(vm.count("SZsymmetry"))
            flatBand->setLambda(vm["SZsymmetry"].as<bool>());
        if(vm.count("hybridizationV"))
            flatBand->setHybridization(vm["hybridizationV"].as<double>());
        if(vm.count("spolarization"))
            flatBand->setSpinPolarization(vm["spolarization"].as<double>());
        chain = std::auto_ptr<nrg::ChainProvider>(flatBand);
    }
    else if(utils::toUpperCase(vm["chain"].as<std::string>()) == "HYBRID") {
        nrg::chain::HybridizationProvider* hybrid = new nrg::chain::HybridizationProvider;
        if(vm.count("config"))
            hybrid->configure(config);
        if(vm.count("PHsymmetry"))
            hybrid->setLambda(vm["PHsymmetry"].as<bool>());
        if(vm.count("SZsymmetry"))
            hybrid->setLambda(vm["SZsymmetry"].as<bool>());
        chain = std::auto_ptr<nrg::ChainProvider>(hybrid);
    }
    else {
        std::cerr << "unknown chain provider: " << vm["chain"].as<std::string>() << std::endl;
        throw std::exception();
    }
    if(vm.count("chainLength"))
        chain->setLength(vm["chainLength"].as<double>());
    if(vm.count("lambda"))
        chain->setLambda(vm["lambda"].as<double>());

    // create & configure broadener
    broadener = std::auto_ptr<nrg::Broadener>(new nrg::broadening::FFTBroadener);
    if(vm.count("config"))
        broadener->configure(config);
    if(vm.count("PPD"))
        broadener->setPolesPerDecade(vm["PPD"].as<int>());
    if(vm.count("peakWidth"))
        broadener->setPeakWidth(vm["peakWidth"].as<double>());

    NRG = std::auto_ptr<nrg::NRG>(new nrg::NRG(*chain.get(), *broadener.get()));
    if(vm.count("config"))
        NRG->configure(config);
    if(vm.count("U"))
        NRG->setU(vm["U"].as<double>());
    if(vm.count("mu"))
        NRG->setEpsF(-NRG->getU()/2.0-vm["mu"].as<double>());
    if(vm.count("T"))
        NRG->setTemperature(vm["T"].as<double>());
    if(vm.count("Ecut"))
        NRG->setEnergyCutOff(vm["Ecut"].as<double>());
    if(vm.count("Ec"))
        NRG->setClusterEnergy(vm["Ec"].as<double>());
    if(vm.count("maxHSdim"))
        NRG->setMaxHilbertSpaceDimension(vm["maxHSdim"].as<int>());



    NRG->init();
    NRG->showInfo();
    NRG->solve();

    std::string outputDir = "";
    if(vm.count("output")) {
        outputDir = vm["output"].as<std::string>();
        if(outputDir[outputDir.size()-1] != '/') {
            outputDir.append("/");
        }
        std::cout << "output directory: " << outputDir << std::endl;
    }

    //output
    if(chain->isSZsymmetric()) {
        math::CFunction S;
        NRG->getSelfEnergy(S);
        S.write(outputDir + "S.txt");
        S.writeBinary(outputDir + "S.dat");
        math::CFunction G;
        NRG->getGreensFunction(G);
        G.write(outputDir + "G.txt");
        G.writeBinary(outputDir + "G.dat");
        math::CFunction F;
        NRG->getFFunction(F);
        F.write(outputDir + "F.txt");
        F.writeBinary(outputDir + "F.dat");
    } else {

    }

}
}
