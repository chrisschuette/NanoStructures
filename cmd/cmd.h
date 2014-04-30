#ifndef CMD_H
#define CMD_H

#include <boost/program_options.hpp>

namespace po = boost::program_options;

po::variables_map parseCommandLine(int argc, char ** argv);


#endif // CMD_H
