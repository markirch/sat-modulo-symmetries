#ifndef SMS_OPTIONS_H
#define SMS_OPTIONS_H

#include "sms.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

extern po::options_description all_opts; // can be used to extend the options

#define DEFAULT_FREQUENCY 30 // default frequency for all propagators if not specified otherwise

typedef struct
{
    bool generate_connected = false;

    bool planar = false;
    int planarityFrequency = DEFAULT_FREQUENCY;

    string forbiddenSubgraphFile;
    string forbiddenInducedSubgraphFile;
    int frequencyForbiddenSubgraphs = DEFAULT_FREQUENCY;

    // coloring related propagators
    int minChromaticNumber = 0;
    int coloringAlgo = 0; // 0 means simple recursive; 1 means simple DPLL-based; 2 means SAT-based
    bool non010colorable = false;
    int triangleVars = 0;

    string qcirFile;
    bool polarityHashing = false;
} propagators_config_t;

/**
 *  Parse the options and create a solver based on the options including predefined propagators.
 *  vm can be use to access the parsed options and modify them.
 *
 */
GraphSolver *parseOptions(const int argc, const char **argv, po::variables_map &vm);

/**
 *  Create a solver based on the configs not necessarily parsed from the command line.
 */
GraphSolver *createSolver(SolverConfig &config, struct minimality_config_t &minimalitConfig, propagators_config_t &propagatorsConfig, string dimacsFile);

#endif