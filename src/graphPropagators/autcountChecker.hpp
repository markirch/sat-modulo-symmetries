#ifndef _AUTCOUNT_CHECKER_H_
#define _AUTCOUNT_CHECKER_H_

#include "../useful.h"
#include "../graphChecker.hpp"
#include "automorphism_counter.hpp"

/**
 * Automorphism counting propagator for SMS
 * Ensures generated graphs have at least k automorphisms
 */
class AutcountChecker : public ComplexFullyDefinedGraphChecker
{
private:
    AutomorphismCounter counter;
    int minAutomorphisms;
    int cutoff;
    bool aggressiveBypass;
    
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable);

public:
    AutcountChecker(int vertices, int minAut = 1, int cutoff = 10000, bool aggressive = false)
        : counter(vertices, false),  // undirected graphs
          minAutomorphisms(minAut),
          cutoff(cutoff),
          aggressiveBypass(aggressive)
    {
        name = "AutcountChecker";
        // ComplexFullyDefinedGraphChecker only runs on final graphs - no frequency logic needed
        (void)cutoff;     // Suppress unused parameter warning
        (void)aggressive; // Suppress unused parameter warning
    }
};

#endif