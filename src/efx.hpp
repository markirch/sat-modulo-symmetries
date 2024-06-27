#ifndef EFX_HPP
#define EFX_HPP
/**
 * Propagator for ensuring acyclicity in the context of EFX (envy free allocation)
 */

#include "useful.h"
#include "graphChecker.hpp"
#include "cadical.hpp"

class EFXPropagator : public PartiallyDefinedGraphChecker
{
private:
    CaDiCaL::Solver solver;             // used for trying to find the cycles
    vector<vector<int>> selectedVertex; // variables indicating the vertices in the cycle
    vector<int> lenghtOfCycle;          // variable indicating the length of the cycle
    vector<vector<int>> edgeVars;       // variables indicating the directed edges in the given graph. Used as assumptions
    int n;                              // number of vertices
    int p;                             // number of partitions
    
    // variables from second version
    vector<int> selectedVertices;
    vector<vector<int>> remainingEdges;
public:
    // p gives the number of partitions (i.e., the number of vertices in the complete directed graph)
    EFXPropagator(int frequency, int n, int p) ;
    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif