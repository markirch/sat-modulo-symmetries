#ifndef FORBIDDEN_SUBGRAPH_H
#define FORBIDDEN_SUBGRAPH_H

#include "graphChecker.hpp"

#ifdef GLASGOW

#include "homomorphism.hh"

class ForbiddenSubgraphCheckerGlasgow : public PartiallyDefinedGraphChecker
{
    vector<InputGraph> forbiddenSubgraphsGlasgow;

    HomomorphismParams defaultHomParams = {
        .timeout = make_shared<Timeout>(0s),
        .restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier),
        .no_supplementals = true};

    InputGraph toGlasgowGraph(const adjacency_matrix_t &adjacencyMatrix, int nVertices)
    {
        InputGraph G(nVertices, false, false);
        for (int i = 0; i < nVertices; i++)
        {
            for (int j = i + 1; j < nVertices; j++)
            {
                if (adjacencyMatrix[i][j] == truth_value_true)
                {
                    G.add_edge(i, j);
                }
            }
        }
        return G;
    }

public:
    ForbiddenSubgraphCheckerGlasgow(int frequency, std::ifstream &forbiddenSubgraphFile)
    {
        checkFinal = false;
        name = "ForbiddenSubgraphCheckerGlasgow";
        this->frequency = frequency;
        // parse file for forbidden subgraphs
        string line;
        while (getline(forbiddenSubgraphFile, line))
        {
            if (strncmp(line.c_str(), "c\t", 2) == 0)
                continue;

            vector<pair<int, int>> edges;
            int max_vertex = 0;

            istringstream iss(line);
            clause_t clause;
            string space_delimiter = " ";

            string lit;
            while (std::getline(iss, lit, ' '))
            {
                int v1 = stoi(lit);
                std::getline(iss, lit, ' ');
                int v2 = stoi(lit);
                edges.push_back(make_pair(v1, v2));
                if (v1 > max_vertex)
                {
                    max_vertex = v1;
                }
                if (v2 > max_vertex)
                {
                    max_vertex = v2;
                }
            }

            forbiddenSubgraphsGlasgow.push_back(InputGraph(max_vertex + 1, false, false));
            for (pair e : edges)
            {
                forbiddenSubgraphsGlasgow.back().add_edge(e.first, e.second);
            }
        }
    }

    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif

#endif // end of header