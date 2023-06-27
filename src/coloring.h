#ifndef COLORING_H
#define COLORING_H

#include "useful.h"

typedef std::vector<int> coloring_t;
typedef std::vector<vertex_t> vertex_ordering_t;

class Coloring
{
private:
    int maxColors;

public:
    Coloring(int maxColors)
    {
        this->maxColors = maxColors;
    }

    bool getColoring(int nVertices, const adjacency_matrix_t &matrix, coloring_t &coloring, int coloringAlgo, vector<int> clique); // get a vertex coloring given an adjacency matrix. Returns true if coloring was found
    static std::vector<edge_t> coloring2monochromaticEdges(const coloring_t &coloring); // get all edges, where vertices have different colors
    static vertex_ordering_t coloring2vertexOrdering(coloring_t &coloring);
};

int get010Coloring(int nVertices, const adjacency_matrix_t &adjacencyMatrix, coloring_t &coloring, vector<vector<vector<int>>> triangle_stats, vector<vector<int>> edge_stats);
int getHyperColoring(int nVertices[2], const adjacency_matrix_t &matrix, vector<int> &coloring);
vector<lit_t> get010ColoringClause(const coloring_t &coloring, int nVertices, const vector<vector<lit_t>> &E, const vector<vector<vector<lit_t>>> &T);
vector<lit_t> getHyperColoringClause(const coloring_t &coloring, int nVertices[2], const adjacency_matrix_t &matrix, const vector<vector<lit_t>> &E);
vector<vector<lit_t>> getHyperColoringCircuit(const coloring_t &coloring, int nVertices[2], const vector<vector<lit_t>> &E, int &next_var);

#endif
