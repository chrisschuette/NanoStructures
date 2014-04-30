#include "cmd.h"

#include <iostream>

po::variables_map parseCommandLine(int argc, char ** argv) {
    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "produce help message")
            ("mpi", "enable message passing interface")
            ("mode", po::value<std::string>(), "nrg / dmft / nano")
            ("config", po::value<std::string>(), "config file")
            ("chain", po::value<std::string>(), "nrg: flat / hybrid")
            ("chainLength", po::value<int>(), "nrg: wilson chain length")
            ("lambda", po::value<double>(), "nrg: logarithmic discretization")
            ("PHsymmetry", po::value<bool>(), "nrg: PH symmetry")
            ("SZsymmetry", po::value<bool>(), "nrg: SZ symmetry")
            ("hybridizationV", po::value<double>(), "flatband: hybridization strength")
            ("spolarization", po::value<double>(), "nrg: spin polarization")
            ("PPD", po::value<int>(), "broadener: poles per decade")
            ("peakWidth", po::value<double>(), "broadener: peak width")
            ("U", po::value<double>(), "nrg/dmft: hubbard U")
            ("mu", po::value<double>(), "nrg/dmft: chemical potential")
            ("T", po::value<double>(), "nrg/dmft/nano: temperature")
            ("Ec", po::value<double>(), "nrg: cluster energy")
            ("Ecut", po::value<double>(), "nrg: energy cut off")
            ("maxHSdim", po::value<int>(), "nrg: max Hilbert space dim")
            ("initialS", po::value<std::string>(), "dmft: initial self energy")
            ("initialD", po::value<std::string>(), "dmft: initial hybridization function")
            ("output", po::value<std::string>(), "nrg/dmft/nano: output directory")
            ("structure", po::value<std::string>(), "nano: structure file")
            ("resume", po::value<std::string>(), "nano: resume directory")
            ("maxIterations", po::value<int>(), "dmft/nano: maximum number of iterations")
            ("doECR", "nano: electronic charge reconstruction")
            ("calcCond", "nano: calculate longitudinal conductivity")
            ("calcHall", "nano: calculate Hall conductivity")
            ("verbose", "nano: verbose")
            ("alphaV", po::value<double>(), "nano: alphaV")
            ("symmetricNano", "nano: PH symmetry under z->Nz-z?")
            ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(desc).run(), vm);
    po::notify(vm);

    // produce help message
    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(0);
    }
    return vm;
}
