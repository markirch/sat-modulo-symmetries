#ifndef SUBGRAPH_ISOMORPHISM_HPP
#define SUBGRAPH_ISOMORPHISM_HPP

#include "useful.h"
#include "graphChecker.hpp"
#include "cadical.hpp"
#include <cstring>

std::pair<cnf_t, vector<vector<int>>> createEncoding(const adjacency_matrix_t &pattern, const int nt);

class ForbiddenSubgraphChecker : public PartiallyDefinedGraphChecker
{
    vector<adjacency_matrix_t> forbiddenSubgraphs; // undefined edges are interpreted as don't cares
    CaDiCaL::Solver forbiddenSubgraphSolver;
    cnf_t forbiddenSubgraphEncoding;
    vector<vector<int>> edgeVars;

public:
    ForbiddenSubgraphChecker(int frequency, std::ifstream &forbiddenSubgraphFile, int n)
    {
        // checkFinal = false;
        name = "ForbiddenSubgraphChecker";
        this->frequency = frequency;
        // parse file for forbidden subgraphs
        if (forbiddenSubgraphFile.is_open())
        {
            load_graphs(forbiddenSubgraphFile);
        }

        if (forbiddenSubgraphs.size() != 1)
            EXIT_UNWANTED_STATE

        auto encoding = createEncoding(forbiddenSubgraphs[0], n);
        forbiddenSubgraphEncoding = encoding.first;
        edgeVars = encoding.second;

        // add encoding to solver
        for (auto &clause : forbiddenSubgraphEncoding)
        {
            for (auto &lit : clause)
            {
                forbiddenSubgraphSolver.add(lit);
            }
            forbiddenSubgraphSolver.add(0);
        }
    }

    void checkProperty(const adjacency_matrix_t &matrix);

    void load_graphs(std::ifstream &graphs, bool induced = false);
};

#endif