#include "useful.h"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "sms.hpp"
#include "options.hpp"

// include all graph propagators


int main(int argc, char const **argv)
{
    clock_t start = clock();

    try
    {
        po::variables_map vm;
        GraphSolver *solver = parseOptions(argc, argv, vm);

        int result = solver->sms_solve();

        delete solver;

        printf("Total time: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
        return result;
    }
    catch (const std::exception &error)
    {
        std::cerr << "Error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
}
