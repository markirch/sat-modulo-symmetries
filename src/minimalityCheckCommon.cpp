#include "minimalityCheck.hpp"

void MinimalityChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    try
    {

        for (auto ordering : vertexOrderings)
        {
            auto matrixCopy = matrix;
            checkMinimality(matrixCopy, ordering, config); // intial partition and vertex ordering are not changed
        }
    }
    catch (LimitReachedException e)
    {
        limit_reached++;
    }
    catch (minimalit_check_result_t e)
    {
        forbidden_graph_t signedEdges = e.clause;
        // TODO currently flipped edges but later replace it in the minimality check itself
        // flip sign of all the edges
        for (size_t i = 0; i < signedEdges.size(); i++)
            if (signedEdges[i].first == truth_value_true)
                signedEdges[i].first = truth_value_false;
            else
                signedEdges[i].first = truth_value_true;

        if (symBreakClauses)
        {
            for (size_t i = 0; i < signedEdges.size(); i++)
                fprintf(symBreakClauses, "%s(%d,%d) ", signedEdges[i].first == truth_value_true ? "-" : "", signedEdges[i].second.first, signedEdges[i].second.second);
            // fprintf(symBreakClauses, "\n");

            fprintf(symBreakClauses, ";");
            // print permutation
            for (auto v: e.permutation)
                fprintf(symBreakClauses, " %d", v);
            fprintf(symBreakClauses, "\n");
        }
        throw signedEdges;
    }
}

void MultipleMinimalityChecker::checkProperty(const vector<adjacency_matrix_t> &matrices)
{
    try
    {
        for (auto ordering : vertexOrderings)
        {
            auto matrixCopy = matrices;
            checkMinimalityMultiple(matrixCopy, ordering, config); // intial partition and vertex ordering are not changed
        }
    }
    catch (LimitReachedException e)
    {
        limit_reached++;
    }
    catch (minimalit_check_result_multi_t e)
    {
        vector<forbidden_graph_t> signedEdgesMulti = e.clause;
        // printf("Graph checked for minimality BEFORE!!:\n");
        // for (int m = 0; m < matrices.size(); m++)
        // {
        //     printf("Layer %d\n", m);
        //     printAdjacencyMatrix(matrices[0], true);
        // }
        // printf("Forbidden:\n");
        // for (int m = 0; m < matrices.size(); m++)
        // {
        //     printf("Layer %d\n", m);
        //     for (size_t i = 0; i < signedEdgesMulti[m].size(); i++)
        //         if (signedEdgesMulti[m][i].first == truth_value_true)
        //             printf("(%d,%d) ", signedEdgesMulti[m][i].second.first, signedEdgesMulti[m][i].second.second);
        //         else
        //             printf("-(%d,%d) ", signedEdgesMulti[m][i].second.first, signedEdgesMulti[m][i].second.second);
        //     printf("\n");
        // }

        // TODO currently flipped edges but later replace it in the minimality check itself
        // flip sign of all the edges
        for (int m = 0; m < (int) matrices.size(); m++) // flip for all of them separately
            for (size_t i = 0; i < signedEdgesMulti[m].size(); i++)
                if (signedEdgesMulti[m][i].first == truth_value_true)
                    signedEdgesMulti[m][i].first = truth_value_false;
                else
                    signedEdgesMulti[m][i].first = truth_value_true;

        if (symBreakClauses)
        {
            EXIT_UNWANTED_STATE // TODO implement
        }
        throw signedEdgesMulti;
    }
}

void MaximalityChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    try
    {

        for (auto ordering : vertexOrderings)
        {
            auto matrixCopy = matrix;
            for (int i = 0; i < (int)matrix.size(); i++)
                for (int j = 0; j < (int)matrix.size(); j++)
                {
                    if (matrix[i][j] == truth_value_true)
                        matrixCopy[i][j] = truth_value_false;
                    if (matrix[i][j] == truth_value_false)
                        matrixCopy[i][j] = truth_value_true;
                }
            checkMinimality(matrixCopy, ordering, config); // intial partition and vertex ordering are not changed
        }
    }
    catch (LimitReachedException e)
    {
        printf("Limit reached\n");
    }
    catch (minimalit_check_result_t e)
    {
        forbidden_graph_t signedEdges = e.clause;
        // TODO currently flipped edges but later replace it in the minimality check itself
        // flip sign of all the edges
        for (size_t i = 0; i < signedEdges.size(); i++)
            if (signedEdges[i].first == truth_value_true)
                signedEdges[i].first = truth_value_false;
            else
                signedEdges[i].first = truth_value_true;

        // flip again because maximized
        for (size_t i = 0; i < signedEdges.size(); i++)
            if (signedEdges[i].first == truth_value_true)
                signedEdges[i].first = truth_value_false;
            else
                signedEdges[i].first = truth_value_true;
        throw signedEdges;
    }
}

void MinimalityCheckerWithStaticPartition::checkProperty(const adjacency_matrix_t &matrix, const vector<truth_value_t> &currentAssignemnt)
{
    try
    {
        int vertices = (int)matrix.size();
        partition_t initialPartition = vector<bool>(vertices, false);
        initialPartition[0] = true;
        for (int v = 1; v < vertices; v++)
            if (currentAssignemnt[staticParitionVars[v - 1][v]] != truth_value_true) // if true then in same partition otherwise seperated
                initialPartition[v] = true;

        vertex_ordering_t ordering;
        for (int i = 0; i < (int)matrix.size(); i++)
            ordering.push_back(i);
        auto matrixCopy = matrix;
        minimalit_check_config_t config;
        config.cutoff = cutoff;
        config.initial_partition = initialPartition;
        checkMinimality(matrixCopy, ordering, config); // intial partition and vertex ordering are not changed
    }
    catch (LimitReachedException e)
    {
        limit_reached++;
    }
    catch (minimalit_check_result_t e)
    {
        vector<clause_t> clauses;
        clause_t clause;
        // TODO currently flipped edges but later replace it in the minimality check itself
        // flip sign of all the edges
        for (size_t i = 0; i < e.clause.size(); i++)
            if (e.clause[i].first == truth_value_true)
            {
                auto edge = e.clause[i].second;
                clause.push_back(-edges[edge.first][edge.second]); // TODO check if plus or minus
            }
            else
            {
                auto edge = e.clause[i].second;
                clause.push_back(edges[edge.first][edge.second]); // TODO check if plus or minus
            }

        for (int v = 0; v < (int)matrix.size(); v++)
        {
            if (e.permutation[v] != v)
                clause.push_back(-staticParitionVars[v][e.permutation[v]]);
        }

        clauses.push_back(clause);
        throw clauses;
    }
}

// TODO swaps with colorings and connected components
/*
checks whether swapping connected components of two colors can decrease the adjacency matrix, unknown edges are assumed to be present for the connected componenets

bool checkSwaps(adjacency_matrix_t &matrix)
{
    int N = intervallsColoring.size(); // last partition is only the single vertex
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
        {
            // get connected components from vertices form color i and j

            vector<int> componenet(vertices, 0); // 0 if not assigned to a component yet

            // get unmarked vertex
            int unmarkedVertex = intervallsColoring[i].first;
            int nComponents = 0;

            do
            {
                nComponents++;

                vector<int> expand;
                componenet[unmarkedVertex] = nComponents;
                expand.push_back(unmarkedVertex);

                while (!expand.empty())
                {
                    int v = expand[expand.size() - 1];
                    expand.pop_back();

                    for (int u = 0; u < vertices; u++)
                    {
                        if (u == v || componenet[u] == nComponents || matrix[v][u] == truth_value_false)
                            continue;

                        if ((u >= intervallsColoring[i].first && u < intervallsColoring[i].second) ||
                            (u >= intervallsColoring[j].first && u < intervallsColoring[j].second))
                        {
                            componenet[u] = nComponents;
                            expand.push_back(u);
                        }
                    }
                }

                // check number of vertices of each color in the connected componenet
                int ni = 0; // vertices with color i
                int nj = 0; // vertices with color j

                for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                    if (componenet[v] == nComponents)
                        ni++;
                for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                    if (componenet[v] == nComponents)
                        nj++;
                // printf("Size of component %d: %d %d\n", c, ni, nj);
                //  if (nj > ni)
                // {
                //     printf("Could be discarded %d %d; %d %d !!!\n", ni, nj, i, j);
                //     return false;
                // }

                if (ni == nj && ni != intervallsColoring[i].second - intervallsColoring[i].first) // check connected components if not all vertices
                {
                    partition_t partition(vertices, false);
                    for (int i = 0; i < N; i++)
                    {
                        partition[intervallsColoring[i].first] = true;
                        partition[intervallsColoring[i].first + 1] = true; // first vertex in own parition.
                    }
                    partition[vertices - 1] = true; // last vertex in own parition

                    // check swapping this two blocks
                    vertex_ordering_t vertex_ordering(vertices);
                    for (int i = 0; i < vertices; i++)
                        vertex_ordering[i] = i;

                    int pos1 = intervallsColoring[i].first;
                    int pos2 = intervallsColoring[j].first;

                    for (int k = 0; k < ni; k++)
                    {
                        // find the next two vertices to swap (in each color seperately)
                        while (componenet[pos1] != nComponents)
                            pos1++;

                        while (componenet[pos2] != nComponents)
                            pos2++;

                        std::swap(vertex_ordering[pos1], vertex_ordering[pos2]);

                        pos1++;
                        pos2++;
                    }
                    try
                    {
                        checkMinimality(matrix, vertex_ordering, {.initial_partition = partition, .cutoff = 0});
                    }
                    catch (std::vector<signed_edge_t> e)
                    {
                        // printf("Not minimal by connected components\n");

                        if (symBreakClausesFile)
                        {
                            fprintf(symBreakClausesFile, "c;%ld;", e.size());

                            // print connected component
                            for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                                if (componenet[v] == nComponents)
                                    fprintf(symBreakClausesFile, "%d ", v);

                            for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                                if (componenet[v] == nComponents)
                                    fprintf(symBreakClausesFile, "%d ", v);
                            fprintf(symBreakClausesFile, ";");

                            // print the two parititions used
                            fprintf(symBreakClausesFile, "%d ; %d;", i, j);
                        }

                        // add information that connected component
                        for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                        {
                            if (componenet[v] != nComponents)
                                continue;
                            for (int u = intervallsColoring[j].first; u < intervallsColoring[j].second; u++)
                            {
                                if (componenet[u] == nComponents)
                                    continue;
                                e.push_back(make_pair(truth_value_true, make_pair(v, u)));
                            }
                        }

                        for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                        {
                            if (componenet[v] == nComponents)
                                continue;
                            for (int u = intervallsColoring[j].first; u < intervallsColoring[j].second; u++)
                            {
                                if (componenet[u] != nComponents)
                                    continue;
                                e.push_back(make_pair(truth_value_true, make_pair(v, u)));
                            }
                        }

                        throw e;
                    }
                }

                // search for unmarked vertex
                unmarkedVertex = -1;
                for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                    if (componenet[v] == 0)
                        unmarkedVertex = v;
                for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                    if (componenet[v] == 0)
                        unmarkedVertex = v;
            } while (unmarkedVertex != -1);
        }

    return true;
}


 */
