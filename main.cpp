#include "err/exception.h"
#include "cmd/cmd.h"
#include "utils/utils.h"
#include "mpi/openmpi.h"

#include "nrg/donrg.h"
#include "dmft/dodmft.h"
#include "nano/donano.h"

#include <iostream>

#include "math/cfunction.h"

using namespace std;

int main(int argc, char** argv)
{
  //try {
        boost::program_options::variables_map vm = parseCommandLine(argc, argv);

        if(vm.count("mode")) {
            std::string mode = utils::toUpperCase(vm["mode"].as<std::string>());
            if(mode == "NRG")
                nrg::doNRG(vm);
            else if(mode == "DMFT")
                dmft::doDMFT(vm);
            else if(mode == "NANO")
                nano::doNano(vm);
            else {
                std::cerr << "unknown mode: " << mode << std::endl;
                throw std::exception();
            }
        } else {
            std::cerr << "mode of operation missing." << std::endl;
            throw std::exception();
        }

  //  } catch(error::Exception & e) {
  //      std::cout << "exception: " << e.what() << std::endl;
  //  } catch(std::exception& e) {
  //      std::cout << "unknown exception: " << e.what() << std::endl;
  //  }
}


/*

    nano::NanoStructure structure(7);
    math::CFunction SM;
    math::CFunction S;
    S.readBinary("final/S_8.dat");
    SM.readBinary("final2/S_20.dat");
    structure.setSelfEnergies(SM);
    structure.setHoppings(1.0,1.0);
    structure.setUs(20.0);
    structure.setU(3, 8.0);
    structure.setU(4, 8.0);
    structure.setU(5, 8.0);
    structure.setSelfEnergy(5, S);
    structure.setSelfEnergy(4, S);
    structure.setSelfEnergy(3, S);
    structure.setMaxIterations(20);
    structure.solve();

    structure.writeBinary(".");

    */

/*
    double U = 20.0;
    math::CFunction test;
    dmft::DMFT dmft;
    dmft.setScaling(100.0);
    //test.createLogGrid(dmft.getScaling(),1.01,1000);
    //for(int i = 0; i < test.getSize(); i++)
    //   test.setValue(i,U/2.0);
    test.readBinary("final/D_20.dat");
    dmft.setHybridizationFunction(test);
    dmft.setStartWithDelta(true);
    dmft.setU(U);
    dmft.setMaxIterations(11);
    dmft.solve();
    dmft.getGreensFunction().writeBinary("final2/G_" + utils::toString(U) + ".dat");
    dmft.getSelfEnergy().writeBinary("final2/S_" + utils::toString(U) + ".dat");
    dmft.getHybridizationFunction().writeBinary("final2/D_" + utils::toString(U) + ".dat");
    dmft.getGreensFunction().write("final2/G_" + utils::toString(U) + ".txt");
    dmft.getSelfEnergy().write("final2/S_" + utils::toString(U) + ".txt");
    dmft.getHybridizationFunction().write("final2/D_" + utils::toString(U) + ".txt");
    return 0;
*/


