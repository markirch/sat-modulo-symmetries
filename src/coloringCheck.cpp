#include <cassert>
#include <numeric>
#include "coloring.h"
#include "coloringCheck.hpp"

void SubgraphChromaticNumberChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    for (int i = 0; i < this->n; i++)
        for (int j = i + 1; j < this->n; j++)
            if (matrix[i][j] == truth_value_unknown)
                return; // not fully defined so can just return

    bool changeInMatrix = previousMatrix.empty(); // if no previous matrix then definetly different.
    for (int i = 0; i < this->n && !changeInMatrix; i++)
        for (int j = i + 1; j < this->n && !changeInMatrix; j++)
            if (matrix[i][j] != previousMatrix[i][j])
            {
                changeInMatrix = true;
                break;
            }

    previousMatrix = matrix;
    Coloring c(minChomaticNumberSubgraph - 1);
    coloring_t coloring(n);
    if (c.getColoring(n, matrix, coloring, coloringAlgo, vector<int>()))
    {
        auto monoChromaticEdges = c.coloring2monochromaticEdges(coloring);
        forbidden_graph_t forbidden_graph;
        for (auto edge : monoChromaticEdges)
            forbidden_graph.push_back(make_pair(truth_value_false, edge));
        throw forbidden_graph;
    }
    // printf("Was not colorable.\n");
}

void MinChromaticNumberChecker::checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &)
{
    Coloring c(chi - 1);
    int nVertices = (int)matrix.size();
    coloring_t coloring(nVertices);
    if (c.getColoring(nVertices, matrix, coloring, coloringAlgo, vector<int>()))
    {
        vector<clause_t> clauses;
        // printf("Colorable\n");
        auto monoChromaticEdges = c.coloring2monochromaticEdges(coloring);
        vector<lit_t> clause = monochromaticEdges2Clause(monoChromaticEdges);
        clauses.push_back(clause);

        if (permuteColorings)
        {
            for (int i = 0; i < (int)coloring.size(); i++)
                for (int j = i + 1; j < (int)coloring.size(); j++)
                {
                    if (coloring[i] == coloring[j])
                        continue;

                    std::swap(coloring[i], coloring[j]); // swap color of vertex i and j
                    auto monoChromaticEdges = c.coloring2monochromaticEdges(coloring);
                    vector<lit_t> clause = monochromaticEdges2Clause(monoChromaticEdges);
                    clauses.push_back(clause);
                    std::swap(coloring[i], coloring[j]);
                }
        }
        throw clauses;
    }
}

void HypergraphMinChromaticNumberChecker::checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &)
{
    Coloring c(chi - 1);
    int nVertices = b_vertices[1]; // size of intersection graph
    coloring_t coloring(nVertices);
    adjacency_matrix_t intersectionMatrix = getIntersectionMatrix(matrix, b_vertices);
    // printf("Size of intersection matrix: %ld\n", intersectionMatrix.size());
    if (c.getColoring(nVertices, intersectionMatrix, coloring, coloringAlgo, getCliqueInIntersectionMatrix(matrix)))
    {
        vector<clause_t> clauses;
        // printf("Colorable\n");
        auto monoChromaticEdges = c.coloring2monochromaticEdges(coloring);
        vector<lit_t> clause = monochromaticEdges2Clause(monoChromaticEdges);
        clauses.push_back(clause);

        if (permuteColorings)
        {
            for (int i = 0; i < (int)coloring.size(); i++)
                for (int j = i + 1; j < (int)coloring.size(); j++)
                {
                    if (coloring[i] == coloring[j])
                        continue;

                    std::swap(coloring[i], coloring[j]); // swap color of vertex i and j
                    auto monoChromaticEdges = c.coloring2monochromaticEdges(coloring);
                    vector<lit_t> clause = monochromaticEdges2Clause(monoChromaticEdges);
                    clauses.push_back(clause);
                    std::swap(coloring[i], coloring[j]);
                }
        }
        throw clauses;
    }
}

vector<int> HypergraphMinChromaticNumberChecker::getCliqueInIntersectionMatrix(const adjacency_matrix_t &incidence_matrix)
{
    adjacency_matrix_t matrix(b_vertices[1], vector<truth_value_t>(b_vertices[1]));
    // compute vertex contained in the most hyperedges
    int maxVertex = 0;
    int maxEdges = 0;
    for (int v = 0; v < b_vertices[0]; v++)
    {
        int counter = 0;
        for (int e = 0; e < b_vertices[1]; e++)
            if (incidence_matrix[v][e + b_vertices[0]] == truth_value_true)
                counter++;
        if (counter > maxEdges)
        {
            maxEdges = counter;
            maxVertex = v;
        }
    }
    vector<int> clique;
    for (int e = 0; e < b_vertices[1]; e++)
        if (incidence_matrix[maxVertex][e + b_vertices[0]] == truth_value_true)
            clique.push_back(e);

    int clique_idx = 0;
    int clique_size = clique.size();
    for (int e = 0; e < b_vertices[1]; e++)
    {
        while (clique_idx < clique_size && clique[clique_idx] < e)
            clique_idx++;
        if (clique_idx < clique_size && e == clique[clique_idx])
        {
            clique_idx++;
            continue;
        }
        bool can_be_added = true;
        for (int f : clique)
        {
            bool intersects = false;
            for (int v = 0; v < b_vertices[0]; v++)
            {
                if (incidence_matrix[v][e + b_vertices[0]] == truth_value_true && incidence_matrix[v][f + b_vertices[0]] == truth_value_true)
                {
                    intersects = true;
                    break;
                }
            }
            if (!intersects)
            {
                can_be_added = false;
                break;
            }
        }
        if (can_be_added)
        {
            clique.push_back(e);
        }
    }

    return clique;
}

#include "cadical.hpp"

void Non010colorableChecker::checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &)
{
    int vertices = (int)matrix.size();
    if (!allColorings)
    {

        coloring_t coloring(vertices);
        if (get010Coloring(vertices, matrix, coloring, *triangle_stats, *edge_stats))
        {
            vector<clause_t> clauses;

            /*for (int v = 0; v < vertices; v++) {
              fprintf(stderr, "%d", coloring[v]);
            }
            fprintf(stderr, "\n");*/
            auto clause = get010ColoringClause(coloring, vertices, edges, triangleVars);
            clauses.push_back(clause);

            if (permuteColorings)
            {
                for (int i = 0; i < vertices; i++)
                    for (int j = i + 1; j < vertices; j++)
                    {
                        if (coloring[i] != coloring[j])
                        {
                            std::swap(coloring[i], coloring[j]);
                            clauses.push_back(get010ColoringClause(coloring, vertices, edges, triangleVars));
                            std::swap(coloring[i], coloring[j]); // undo
                        }
                    }
            }
            throw clauses;
        }
    }
    else // get all valid colorings
    {
        auto coloringSolver = new CaDiCaL::Solver();
        vector<int> color(vertices);
        iota(color.begin(), color.end(), 1);

        for (int i = 0; i < vertices; i++)
            for (int j = i + 1; j < vertices; j++)
                if (matrix[i][j] == truth_value_true)
                {
                    coloringSolver->add(color[i]);
                    coloringSolver->add(color[j]);
                    coloringSolver->add(0);
                }

        for (int i = 0; i < vertices; i++)
            for (int j = i + 1; j < vertices; j++)
                for (int k = j + 1; k < vertices; k++)
                    if (matrix[i][j] == truth_value_true && matrix[i][k] == truth_value_true && matrix[k][j] == truth_value_true)
                    {
                        coloringSolver->add(-color[i]);
                        coloringSolver->add(-color[j]);
                        coloringSolver->add(-color[k]);
                        coloringSolver->add(0);
                    }

        int numberOfColorings = 0;
        vector<clause_t> clauses;

        while (coloringSolver->solve() == 10)
        {
            numberOfColorings++;
            vector<int> resultingColoring(vertices, 0);
            for (int i = 0; i < vertices; i++)
                if (coloringSolver->val(color[i]) < 0)
                    resultingColoring[i] = 1;

            clauses.push_back(get010ColoringClause(resultingColoring, vertices, edges, triangleVars));

            for (int i = 0; i < vertices; i++)
                if (resultingColoring[i] == 1)
                    coloringSolver->add(color[i]);
            coloringSolver->add(0);
        }
        printf("Number of colorings: %d\n", numberOfColorings);
        delete coloringSolver;

        if (numberOfColorings != 0)
            throw clauses;
    }
}

void HyperColoringChecker::checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &nextFreeVariable)
{
    coloring_t coloring(b_vertices[0]);
    if (getHyperColoring(b_vertices, matrix, coloring))
    {
        throw getHyperColoringCircuit(coloring, b_vertices, edges, nextFreeVariable);
    }
}

void GreedyColoring::checkProperty(const adjacency_matrix_t &matrix, const vector<int> &, int &nextFreeLit)
{
    Coloring c(chi - 1);
    int nVertices = (int)matrix.size();
    coloring_t coloring(nVertices);
    // printf("chromatic number chi: %d\n", chi);
    if (c.getColoring(nVertices, matrix, coloring, coloringAlgo, vector<int>()))
    {
        // printf("Found chi - 1 coloring\n");
        vector<clause_t> clauses;
        // printf("Colorable\n");
        // TODO create clauses with greedy ordering

        // printf("START\n");
        int n = (int)matrix.size();
        vector<int> vertexOrdering(n);

        int pos = 0;
        for (int c = 0; c < chi - 1; c++)
            for (int i = 0; i < n; i++)
                if (coloring[i] == c)
                    vertexOrdering[pos++] = i;

        assert(pos == n);

        // coloring[i][0] denotes if vertex vertexOrdering[i] has color 0
        vector<vector<lit_t>> coloring(n, vector<int>(chi - 1));  // coloring of each vertex
        vector<vector<lit_t>> available(n, vector<int>(chi - 1)); // check if color is available, i.e., no smaller vertex in ordering has the color
        vector<lit_t> isColored(n);          // check if vertex is colored
        vector<vector<vector<lit_t>>> adjacentAndColorC(n, vector<vector<int>>(n, vector<int>(chi-1)));

        for (int i = 0; i < n; i++)
        {
            isColored[i] = nextFreeLit++;
            //printf("isColored %d  %d\n", i, isColored[i]);
            for (int c = 0; c < chi - 1; c++)
            {
                coloring[i][c] = nextFreeLit++;
                // printf("Coloring %d %d: %d\n", i, c, coloring[i][c]);
                available[i][c] = nextFreeLit++;
                // printf("available %d %d: %d\n", i, c, available[i][c]);
                for (int j = i + 1; j < n; j++)
                {
                    adjacentAndColorC[i][j][c] = nextFreeLit++;
                    // printf("adjacentAndColorC[%d][%d][%d]: %d\n", i, j, c, adjacentAndColorC[i][j][c]);
                }
            }
        }

        // --------------add clauses-----------

        // one color implies isColored
        for (int i = 0; i < n; i++)
            for (int c = 0; c < chi - 1; c++)
                clauses.push_back({-coloring[i][c], isColored[i]});

        // at most on color
        for (int i = 0; i < n; i++)
            for (int c1 = 0; c1 < chi - 1; c1++)
                for (int c2 = c1 + 1; c2 < chi - 1; c2++)
                {
                    clauses.push_back({-coloring[i][c1], -coloring[i][c2]});
                }

        // at least one node is not colored
        vector<int> clause(n);
        for (int i = 0; i < n; i++)
            clause[i] = -isColored[i];
        clauses.push_back(clause);

        // smallest available means that it should get this color
        for (int i = 0; i < n; i++)
        {
            for (int c1 = 0; c1 < chi - 1; c1++)
            {
                vector<int> clause;
                for (int c2 = 0; c2 < c1; c2++) // smaller colors
                    clause.push_back(available[i][c2]);

                clause.push_back(-available[i][c1]);
                clause.push_back(coloring[i][c1]);
                clauses.push_back(clause);
            }
        }

        // truth_value_true and color c
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                for (int c = 0; c < chi - 1; c++)
                {

                    int v1 = std::min(vertexOrdering[i], vertexOrdering[j]);
                    int v2 = std::max(vertexOrdering[i], vertexOrdering[j]);
                    clauses.push_back({edges[v1][v2], -adjacentAndColorC[i][j][c]});  // not truth_value_true implies not truth_value_true and color c
                    clauses.push_back({coloring[i][c], -adjacentAndColorC[i][j][c]}); // not truth_value_true implies not truth_value_true and color c
                    clauses.push_back({-edges[v1][v2], -coloring[i][c], adjacentAndColorC[i][j][c]});
                }
            }
        }

        // if no smaller vertex is truth_value_true and has color c, than color c is available.
        for (int i = 0; i < n; i++)
        {
            for (int c = 0; c < chi - 1; c++)
            {
                vector<lit_t> clause;
                for (int j = 0; j < i; j++)
                    clause.push_back(adjacentAndColorC[j][i][c]);
                clause.push_back(available[i][c]);
                clauses.push_back(clause);
            }
        }

        for (int i = 0; i < n; i++)
            for (int c = 0; c < chi - 1; c++)
                for (int j = 0; j < i; j++)
                {
                    // TODO maybe use edge literals and color literals directly.
                    clauses.push_back({-adjacentAndColorC[j][i][c], -available[i][c]}); // if truth_value_true and color c than not available
                }

        // clauses.push_back({});

        // for (auto vec : clauses)
        // {
        //     for (const auto &element : vec)
        //     {
        //         std::cout << element << " ";
        //     }
        //     std::cout << std::endl;
        // }

        throw clauses;
    }
}
