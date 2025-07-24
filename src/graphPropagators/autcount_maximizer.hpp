#ifndef _AUTCOUNT_MAXIMIZER_H_
#define _AUTCOUNT_MAXIMIZER_H_

#include "../useful.h"
#include "../graphChecker.hpp"
#include "automorphism_counter.hpp"

// Forward declaration for GraphSolver
class GraphSolver;

/**
 * Adaptive automorphism maximizer for SMS (final graphs only)
 * Tracks and updates the maximum automorphism count found during search
 */
class AutcountMaximizerChecker : public ComplexFullyDefinedGraphChecker
{
private:
    AutomorphismCounter counter;
    GraphSolver& solver;
    int cutoff;
    
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable);

public:
    AutcountMaximizerChecker(GraphSolver& s, int vertices, int cutoff_val = 10000)
        : solver(s), counter(vertices, false), cutoff(cutoff_val)
    {
        name = "AutcountMaximizerChecker";
    }
};

/**
 * Adaptive automorphism maximizer for SMS (partial graphs during search)
 * Uses dynamic threshold from solver state to prune search space
 */
class AutcountMaximizerPDGChecker : public PartiallyDefinedGraphChecker
{
private:
    AutomorphismCounter counter;
    GraphSolver& solver;
    int maxUndefinedEdges;
    
    void checkProperty(const adjacency_matrix_t &matrix);

public:
    AutcountMaximizerPDGChecker(GraphSolver& s, int vertices, int frequency = DEFAULT_FREQ, int maxUndefined = 15)
        : PartiallyDefinedGraphChecker(frequency),
          solver(s), 
          counter(vertices, false),
          maxUndefinedEdges(maxUndefined)
    {
        name = "AutcountMaximizerPDGChecker";
    }
};

#endif