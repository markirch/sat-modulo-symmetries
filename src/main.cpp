#include "useful.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "sms.hpp"
#include "options.hpp"

// include all graph propagators


int main(int argc, char const **argv)
{
    clock_t start = clock();

    po::variables_map vm;
    GraphSolver *solver = parseOptions(argc, argv, vm);

    // TODO parse arguments
    solver->sms_solve();

    delete solver;

    printf("Total time: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
    return 0;
}
