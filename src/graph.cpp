#include "graph.hpp"

void printAdjacencyMatrix(const adjacency_matrix_t &matrix, const bool directed, const bool printFullMatrix)
{
    if (printFullMatrix)
    {
        for (auto row : matrix)
        {
            for (auto value : row)
                printf("%d ", value);
            printf("\n");
        }
    }
    else
    {
        printf("[");
        bool first = true;

        for (int i = 0; i < (int)matrix.size(); i++)
        {
            for (int j = directed ? 0 : i + 1; j < (int)matrix.size(); j++)
            {
                if (directed && i == j)
                    continue;
                if (matrix[i][j] == truth_value_true)
                {
                    if (!first)
                        printf(",");
                    printf("(%d,%d)", i, j);
                    first = false;
                }
            }
        }
        printf("]\n");
    }
}

GraphHandler::GraphHandler(int vertices, bool directed) : vertices(vertices), directed(directed)
{
    if (directed)
        numVariables = vertices * vertices - vertices;
    else
        numVariables = vertices * (vertices - 1) / 2;

    if (overlayingGraphs > 1)
    {
        numVariables *= overlayingGraphs;
        EXIT_UNWANTED_STATE // not supported yet
    }

    vertexPair2Lit = vector<vector<lit_t>>(vertices, vector<lit_t>(vertices, 0));
    int counter = 1;
    for (int i = 0; i < vertices; i++)
    {
        for (int j = directed ? 0 : i + 1; j < vertices; j++)
        {
            if (directed && i == j)
                continue;
            vertexPair2Lit[i][j] = counter++;

            if (!directed)
                vertexPair2Lit[j][i] = vertexPair2Lit[i][j];
        }
    }
}

vector<lit_t> GraphHandler::graph2assignment(const adjacency_matrix_t &matrix)
{
    vector<lit_t> assignment;
    for (int i = 0; i < vertices; i++)
    {
        for (int j = directed ? 0 : i + 1; j < vertices; j++)
        {
            if (matrix[i][j] == truth_value_true)
                assignment.push_back(vertexPair2Lit[i][j]);
            else if (matrix[i][j] == truth_value_false)
                assignment.push_back(-vertexPair2Lit[i][j]);
        }
    }
    return assignment;
}

adjacency_matrix_t GraphHandler::assignment2graph(const vector<truth_value_t> &assignment)
{
    adjacency_matrix_t matrix(vertices, vector<truth_value_t>(vertices, truth_value_unknown));
    for (int i = 0; i < vertices; i++)
    {
        for (int j = directed ? 0 : i + 1; j < vertices; j++)
        {
            if (assignment[vertexPair2Lit[i][j]] == truth_value_true)
                matrix[i][j] = truth_value_true;
            else if (assignment[vertexPair2Lit[i][j]] == truth_value_false)
                matrix[i][j] = truth_value_false;

            if (!directed)
            {
                // make symmetric
                matrix[j][i] = matrix[i][j];
            }
        }
    }
    return matrix;
}

clause_t GraphHandler::theClauseThatBlocks(const forbidden_graph_t &fg)
{
    clause_t clause;
    for (auto signedEdge : fg)
    {
        auto edge = signedEdge.second;
        if (signedEdge.first == truth_value_true)
            clause.push_back(-vertexPair2Lit[edge.first][edge.second]);
        else // assum that not truth_value_unknown
            clause.push_back(vertexPair2Lit[edge.first][edge.second]);
    }
    return clause;
}

clause_t GraphHandler::solutionBlockingClause(const vector<truth_value_t> &assignment)
{
    clause_t clause;
    for (int i = 1; i <= numVariables; i++)
    {
        if (assignment[i] == truth_value_true)
            clause.push_back(-i);
        else if (assignment[i] == truth_value_false)
            clause.push_back(i);
    }
   
    return clause;
}

int GraphHandler::numAssigned(const vector<truth_value_t> &assignment)
{
    int num = 0;
    for (int i = 1; i <= numVariables; i++)
    {
        if (assignment[i] != truth_value_unknown)
            num++;
    }
    return num;
}

void GraphHandler::printCube(const vector<truth_value_t> &assignment)
{
    printf("a");
    for (int i = 1; i <= numVariables; i++)
    {
        if (assignment[i] == truth_value_true)
            printf(" %d", i);
        else if (assignment[i] == truth_value_false)
            printf(" %d", -i);
    }
    printf(" 0\n");
}