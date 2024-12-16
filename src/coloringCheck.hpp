#ifndef COLORING_CHECKER_H
#define COLORING_CHECKER_H

#include "graphChecker.hpp"
#include "useful.h"

/**
 * @brief Check wether the subgraph given by the first few vertices has at least a certain chromatic number
 *
 */
class SubgraphChromaticNumberChecker : public PartiallyDefinedGraphChecker
{
    int n; // size of subgraph to check
    int minChomaticNumberSubgraph;
    adjacency_matrix_t previousMatrix;
    int coloringAlgo;
    bool permuteColorings = false;

public:
    SubgraphChromaticNumberChecker(int frequency, int sizeOfSubgraph, int minChomaticNumberSubgraph, int coloringAlgo)
    {
        this->name = "SubgraphChromaticNumberChecker";
        this->frequency = frequency;
        this->n = sizeOfSubgraph;
        this->minChomaticNumberSubgraph = minChomaticNumberSubgraph;
        this->coloringAlgo = coloringAlgo;
    }
    void checkProperty(const adjacency_matrix_t &matrix);
};

class MinChromaticNumberChecker : public ComplexFullyDefinedGraphChecker
{
    int chi; // minimal chromatic number
    int coloringAlgo;
    bool permuteColorings = false;
    vector<vector<int>> edges; // edge variables

    clause_t monochromaticEdges2Clause(vector<edge_t> mches)
    {
        vector<lit_t> clause;
        for (auto e : mches)
            clause.push_back(edges[e.first][e.second]);
        return clause;
    }

public:
    MinChromaticNumberChecker(int minChomaticNumber, int coloringAlgo, vector<vector<int>> edges, bool permuteColorings)
    {
        this->name = "MinChromaticNumberChecker";
        chi = minChomaticNumber;
        this->coloringAlgo = coloringAlgo;
        this->edges = edges;
        addsOnlyObservedLiterals = true;
        this->permuteColorings = permuteColorings;
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &);
};

class HypergraphMinChromaticNumberChecker : public ComplexFullyDefinedGraphChecker
{
    int chi; // minimal chromatic number
    int coloringAlgo;
    bool permuteColorings = false;
    vector<vector<int>> edges_intersection_graph; // variables representing the intersection graph
    int b_vertices[2];                            // size of the vertices of the bipartite graph representing the incidence matrix

    clause_t monochromaticEdges2Clause(vector<edge_t> mches)
    {
        vector<lit_t> clause;
        for (auto e : mches)
        {
            clause.push_back(edges_intersection_graph[e.first][e.second]);
        }
        return clause;
    }

    // given the incidence matrix, the function return a clique of the intersection matrix
    vector<int> getCliqueInIntersectionMatrix(const adjacency_matrix_t &incidence_matrix);

public:
    HypergraphMinChromaticNumberChecker(int minChomaticNumber, int coloringAlgo, vector<vector<int>> edges_intersection_graph, int b_vertices[2])
    {
        this->name = "HypergraphMinChromaticNumberChecker";
        chi = minChomaticNumber;
        this->coloringAlgo = coloringAlgo;
        addsOnlyObservedLiterals = false; // intersection graph is not observed
        this->edges_intersection_graph = edges_intersection_graph;
        this->b_vertices[0] = b_vertices[0];
        this->b_vertices[1] = b_vertices[1];
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &);
};

class Non010colorableChecker : public ComplexFullyDefinedGraphChecker
{
    vector<vector<vector<int>>>& triangleVars; // the variables indicating whether certain triangles are present.
    vector<vector<int>>& edges;
    bool allColorings = false;
    bool permuteColorings = false;

    vector<vector<vector<int>>> *triangle_stats;
    vector<vector<int>> *edge_stats;

public:
    Non010colorableChecker(vector<vector<vector<int>>> &triangleVars, vector<vector<int>> &edges, vector<vector<vector<int>>> *triangle_stats, vector<vector<int>> *edge_stats) :
      triangleVars (triangleVars),
      edges (edges)
    {
        this->triangle_stats = triangle_stats;
        this->edge_stats = edge_stats;
        this->name = "Non010colorableChecker";
        addsOnlyObservedLiterals = false;
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &);
};

// TODO not sure what exactly it does
class HyperColoringChecker : public ComplexFullyDefinedGraphChecker
{
    int coloringAlgo;
    int b_vertices[2]; // size of the vertices of the bipartite graph representing the incidence matrix
    edge_vars_t edges;

public:
    HyperColoringChecker(int coloringAlgo, int b_vertices[2], edge_vars_t edges)
    {
        this->name = "HyperColoringChecker";
        this->coloringAlgo = coloringAlgo;
        addsOnlyObservedLiterals = false;
        this->b_vertices[0] = b_vertices[0];
        this->b_vertices[1] = b_vertices[1];
        this->edges = edges;
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &);
};

class GreedyColoring : public ComplexFullyDefinedGraphChecker
{
    int coloringAlgo;
    edge_vars_t edges;
    int chi; // min chromatic number

public:
    GreedyColoring(int coloringAlgo, edge_vars_t edges, int chi)
    {
        this->name = "Greedy colorings";
        this->coloringAlgo = coloringAlgo;
        addsOnlyObservedLiterals = false;
        this->edges = edges;
        this->chi = chi;
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &);
};
#endif
