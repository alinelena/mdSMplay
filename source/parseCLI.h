
#ifndef PARSECLI_H
#define PARSECLI_H
#include <boost/program_options.hpp>
#include <string>
#include "generalinputs.h"
boost::program_options::variables_map parseCLiandConfig(int ac, char** av, generalInputs *inps);
void printOptions(boost::program_options::variables_map *clioMap, std::string *dbgf);
#endif //PARSECLI_H
