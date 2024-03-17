#ifndef FORBIDDEN_SUBGRAPH_H
#define FORBIDDEN_SUBGRAPH_H

#include "graphChecker.hpp"

#ifdef GLASGOW

#include "gss/homomorphism.hh"
#include <chrono>
using namespace std::chrono;
using namespace gss;
using std::pair;

InputGraph toGlasgowGraph(const adjacency_matrix_t &adjacencyMatrix, int nVertices);
InputGraph toLabeledGlasgowGraph(const adjacency_matrix_t &adjacencyMatrix, int nVertices);

void HomomorphismResultToSubgraph(HomomorphismResult &res, forbidden_graph_t &forbidden_graph, const InputGraph &H);
void InducedHomomorphismResultToSubgraph(HomomorphismResult &res, forbidden_graph_t &forbidden_graph, const InputGraph &H);

extern HomomorphismParams defaultHomParams;

#define PRESENT_LABEL "1"
#define ABSENT_LABEL  "0"
#define UNKNOWN_LABEL  "u"

class MaxCliqueChecker : public PartiallyDefinedGraphChecker
{
    int maxCliqueSize = 0;

public:
    MaxCliqueChecker(int frequency, int maxClique)
    {
        name = "MaxCliqueChecker";
        this->frequency = frequency;
        this->maxCliqueSize = maxClique;
    }

    void checkProperty(const adjacency_matrix_t &matrix)
    {
        InputGraph clique = InputGraph(maxCliqueSize + 1, false, false); // forbids all cliques larger than maxCliqueSize
        for (int i = 0; i < maxCliqueSize + 1; i++)
            for (int j = i + 1; j < maxCliqueSize + 1; j++)
                clique.add_edge(i, j);

        InputGraph G = toGlasgowGraph(matrix, (int)matrix.size());
        HomomorphismResult res = solve_homomorphism_problem(clique, G, defaultHomParams);
        if (!res.mapping.empty())
        {
            forbidden_graph_t forbidden_graph;
            HomomorphismResultToSubgraph(res, forbidden_graph, clique);
            throw forbidden_graph;
        }
    }
};

class MaxIndependentSetChecker : public PartiallyDefinedGraphChecker
{
    int maxIndependentSetSize = 0;

public:
    MaxIndependentSetChecker(int frequency, int maxIndependentSetSize)
    {
        name = "MaxIndependentSetChecker";
        this->frequency = frequency;
        this->maxIndependentSetSize = maxIndependentSetSize;
    }

    void checkProperty(const adjacency_matrix_t &matrix)
    {
        //  negate matrix
        adjacency_matrix_t negatedMatrix = matrix;
        for (int i = 0; i < (int)matrix.size(); i++)
            for (int j = 0; j < (int)matrix.size(); j++)
                if (matrix[i][j] == truth_value_true)
                    negatedMatrix[i][j] = truth_value_false;
                else if (matrix[i][j] == truth_value_false)
                    negatedMatrix[i][j] = truth_value_true;

        InputGraph indSet = InputGraph(this->maxIndependentSetSize + 1, false, false); // forbids all independent set larger than maxIndependentSetSize
        for (int i = 0; i < maxIndependentSetSize + 1; i++)
            for (int j = i + 1; j < maxIndependentSetSize + 1; j++)
                indSet.add_edge(i, j);

        InputGraph G = toGlasgowGraph(negatedMatrix, (int)negatedMatrix.size());
        HomomorphismResult res = solve_homomorphism_problem(indSet, G, defaultHomParams);
        if (!res.mapping.empty())
        {
            forbidden_graph_t forbidden_graph;
            HomomorphismResultToSubgraph(res, forbidden_graph, indSet);
            // negate forbidden graph
            for (int i = 0; i < (int)forbidden_graph.size(); i++)
                if (forbidden_graph[i].first == truth_value_true)
                    forbidden_graph[i].first = truth_value_false;
                else if (forbidden_graph[i].first == truth_value_false)
                    forbidden_graph[i].first = truth_value_true;
            throw forbidden_graph;
        }
    }
};

class ForbiddenSubgraphCheckerGlasgow : public PartiallyDefinedGraphChecker
{
    vector<InputGraph> forbiddenSubgraphsGlasgow;
    vector<InputGraph> forbiddenInducedSubgraphsGlasgow;

    HomomorphismParams inducedHomParams = {
        .timeout = std::make_shared<Timeout>(0s),
        .induced = true,
        .restarts_schedule = std::make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier),
        .no_supplementals = true};

public:
    ForbiddenSubgraphCheckerGlasgow(int frequency, std::ifstream &forbiddenSubgraphFile, std::ifstream &forbiddenInducedSubgraphFile)
    {
        // checkFinal = false;
        name = "ForbiddenSubgraphCheckerGlasgow";
        this->frequency = frequency;
        // parse file for forbidden subgraphs
        if (forbiddenSubgraphFile.is_open())
        {
            load_graphs(forbiddenSubgraphsGlasgow, forbiddenSubgraphFile);
        }
        if (forbiddenInducedSubgraphFile.is_open())
        {
            load_graphs(forbiddenInducedSubgraphsGlasgow, forbiddenInducedSubgraphFile, true);
        }
    }

    void load_graphs(std::vector<InputGraph> &storage, std::ifstream &graphs, bool induced = false)
    {
        string line;
        while (getline(graphs, line))
        {
            if (strncmp(line.c_str(), "c\t", 2) == 0)
                continue;

            vector<pair<int, int>> edges;

            std::istringstream iss(line);
            clause_t clause;
            string space_delimiter = " ";

            string lit;
            std::getline(iss, lit, ' ');
            int n = stoi(lit); // first integer gives number of vertices the remaining the edges

            while (std::getline(iss, lit, ' '))
            {
                int v1 = stoi(lit);
                std::getline(iss, lit, ' ');
                int v2 = stoi(lit);
                edges.push_back(std::make_pair(v1, v2));
            }

            storage.push_back(InputGraph(n, false, induced)); // in the induced case we have edge labels
            if (induced)
            {
                adjacency_matrix_t adjacencyMatrix(n, vector<truth_value_t>(n, truth_value_false));
                for (pair e : edges)
                    adjacencyMatrix[e.second][e.first] = adjacencyMatrix[e.first][e.second] = truth_value_true;

                for (int i = 0; i < n; i++)
                    for (int j = i + 1; j < n; j++)
                    {
                        if (adjacencyMatrix[i][j] == truth_value_true)
                        {
                            storage.back().add_directed_edge(i, j, PRESENT_LABEL);
                            storage.back().add_directed_edge(j, i, PRESENT_LABEL);
                        }
                        else
                        {
                            storage.back().add_directed_edge(i, j, ABSENT_LABEL);
                            storage.back().add_directed_edge(j, i, ABSENT_LABEL);
                        }
                    }
            }
            else
                for (pair e : edges)
                    storage.back().add_edge(e.first, e.second);
        }
    }

    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif // end of GLASGOW

#endif // end of header
