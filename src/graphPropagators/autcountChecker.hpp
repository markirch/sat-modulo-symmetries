#ifndef _AUTCOUNT_CHECKER_H_
#define _AUTCOUNT_CHECKER_H_

#include "../useful.h"
#include "../graphChecker.hpp"
#include "automorphism_counter.hpp"

/**
 * Automorphism counting propagator for SMS (final graphs only)
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

/**
 * Automorphism counting propagator for SMS (partial graphs during search)
 * Prunes search space when PDGs cannot achieve minimum automorphisms
 */
class AutcountPDGChecker : public PartiallyDefinedGraphChecker
{
private:
    AutomorphismCounter counter;
    int minAutomorphisms;
    
    void checkProperty(const adjacency_matrix_t &matrix);

public:
    AutcountPDGChecker(int vertices, int minAut = 1, int frequency = DEFAULT_FREQ)
        : PartiallyDefinedGraphChecker(frequency),
          counter(vertices, false),  // undirected graphs
          minAutomorphisms(minAut)
    {
        name = "AutcountPDGChecker";
        // This checker runs during search on partial graphs
    }
};

#endif