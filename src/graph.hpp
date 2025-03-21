/**
 * Contains the information how the graph is represented and how to convert between the graph and the SAT problem.
 */

#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "useful.h"

typedef std::vector<std::vector<truth_value_t>> adjacency_matrix_t;

// Print the adjacency matrix of a graph. If printFullMatrix is true, print the full matrix, otherwise only the edges as a python list.
void printAdjacencyMatrix(const adjacency_matrix_t &matrix, const bool directed, const bool printFullMatrix);

enum GraphFormat
{
    python_edge_list,
    full_matrix,
    graph6
};

// TODO maybe find better name
class GraphHandler
{
private:
    int vertices;
    int numVariables; // the number of variables to encode the combinatorial object (assuming that the first variables are used for encoding)
    bool directed;
    int overlayingGraphs = 1; // number of overlaying graphs
    enum GraphFormat graphFormat = python_edge_list;

    vector<vector<lit_t>> vertexPair2Lit; // for each edge (i,j) the corresponding variable

public:
    GraphHandler(int vertices, bool directed);

    vector<lit_t> graph2assignment(const adjacency_matrix_t &matrix);             // Given a (partially defined) graph, return the corresponding literals of the defined edges
    adjacency_matrix_t assignment2graph(const vector<truth_value_t> &assignment); // Given an assignment, return the corresponding partially defined graph
    clause_t theClauseThatBlocks(const forbidden_graph_t &fg);                    // Given a forbidden graph, return the clause that blocks it
    clause_t solutionBlockingClause(const vector<truth_value_t> &assignment);     // Given an assignment, return the clause that blocks the corresponding (partially defined) graph

    void print(const adjacency_matrix_t matrix)
    {
        switch (graphFormat)
        {
        case python_edge_list:
            printAdjacencyMatrix(matrix, directed, false);
            break;
        case full_matrix:
            printAdjacencyMatrix(matrix, directed, true);
            break;
        case graph6:
            printf("ERROR: graph 6 format not implemented yet");
            EXIT_UNWANTED_STATE
        default:
            EXIT_UNWANTED_STATE
        }
    }

    int getNumVariables() { return numVariables; }
    vector<vector<lit_t>> getEdgeVariables() { return vertexPair2Lit; }
    int numAssigned(const vector<truth_value_t> &assignment); // number of assigned edge variables
    void printCube(const vector<truth_value_t> &assignment);
};

#endif