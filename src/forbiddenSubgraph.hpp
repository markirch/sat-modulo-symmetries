#ifndef FORBIDDEN_SUBGRAPH_H
#define FORBIDDEN_SUBGRAPH_H

#include "graphChecker.hpp"

#ifdef GLASGOW

#include "gss/homomorphism.hh"
#include <chrono>
using namespace std::chrono;
using namespace gss;
using std::pair;

class ForbiddenSubgraphCheckerGlasgow : public PartiallyDefinedGraphChecker
{
    vector<InputGraph> forbiddenSubgraphsGlasgow;
    vector<InputGraph> forbiddenInducedSubgraphsGlasgow;

    HomomorphismParams defaultHomParams = {
        .timeout = std::make_shared<Timeout>(0s),
        .restarts_schedule = std::make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier),
        .no_supplementals = true
    };

    HomomorphismParams inducedHomParams = {
        .timeout = std::make_shared<Timeout>(0s),
        .induced = true,
        .restarts_schedule = std::make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier),
        .no_supplementals = true
    };

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
    ForbiddenSubgraphCheckerGlasgow(int frequency, std::ifstream &forbiddenSubgraphFile, std::ifstream &forbiddenInducedSubgraphFile)
    {
        //checkFinal = false;
        name = "ForbiddenSubgraphCheckerGlasgow";
        this->frequency = frequency;
        // parse file for forbidden subgraphs
        if (forbiddenSubgraphFile.is_open()) {
          load_graphs(forbiddenSubgraphsGlasgow, forbiddenSubgraphFile);
        }
        if (forbiddenInducedSubgraphFile.is_open()) {
          load_graphs(forbiddenInducedSubgraphsGlasgow, forbiddenInducedSubgraphFile);
        }
    }

    void load_graphs(std::vector<InputGraph> &storage, std::ifstream &graphs) {
        string line;
        while (getline(graphs, line))
        {
            if (strncmp(line.c_str(), "c\t", 2) == 0)
                continue;

            vector<pair<int, int>> edges;
            int max_vertex = 0;

            std::istringstream iss(line);
            clause_t clause;
            string space_delimiter = " ";

            string lit;
            while (std::getline(iss, lit, ' '))
            {
                int v1 = stoi(lit);
                std::getline(iss, lit, ' ');
                int v2 = stoi(lit);
                edges.push_back(std::make_pair(v1, v2));
                if (v1 > max_vertex)
                {
                    max_vertex = v1;
                }
                if (v2 > max_vertex)
                {
                    max_vertex = v2;
                }
            }

            storage.push_back(InputGraph(max_vertex + 1, false, false));
            for (pair e : edges)
            {
                storage.back().add_edge(e.first, e.second);
            }
        }
    }

    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif

#endif // end of header
