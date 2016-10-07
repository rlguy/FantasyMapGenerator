#ifndef GEN_CONFIG_H
#define GEN_CONFIG_H

#include <stdio.h>
#include <iostream>
#include <sstream>

#include "argtable3/argtable3.h"

namespace gen {
namespace config {

struct OptionArgs {
	struct arg_lit *help;
	struct arg_lit *timeseed;
	struct arg_int *seed;
	struct arg_dbl *resolution;
	struct arg_file *outfile;
	struct arg_file *output;
    struct arg_dbl *eroamount;
    struct arg_int *erosteps;
	struct arg_lit *drawinfo;
    struct arg_lit *verbose;
	struct arg_end *end;
};

extern unsigned int seed;
extern double resolution;
extern std::string outfileExt;
extern std::string outfile;
extern double erosionAmount;
extern int erosionIterations;
extern bool verbose;

template<class T>
std::string toString(T item) {
    std::ostringstream sstream;
    sstream << item;

    return sstream.str();
}

void print(std::string msg);
bool parseOptions(int argc, char **argv);
bool _displayInfo(OptionArgs opts);
bool _displayDrawSupportInfo(arg_lit *drawinfo);
bool _setOptions(OptionArgs opts);
bool _setSeed(arg_lit *timeseed, arg_int *seed);
bool _setResolution(arg_dbl *res);
bool _setOutputFile(arg_file *outfile1, arg_file *outfile2);
bool _setErosionAmount(arg_dbl *amount);
bool _setVerbosity(arg_lit *verbose);
bool _setErosionIterations(arg_int *iterations);

}
}

#endif