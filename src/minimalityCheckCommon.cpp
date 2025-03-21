#include "minimalityCheck.hpp"

inline adjacency_matrix_t flipMatrix(const adjacency_matrix_t &matrix)
{
    adjacency_matrix_t m = matrix; // copy

    for (int i = 0; i < (int)m.size(); i++)
        for (int j = 0; j < (int)m.size(); j++)
        {
            if (matrix[i][j] == truth_value_true)
                m[i][j] = truth_value_false;
            else if (matrix[i][j] == truth_value_false)
            {
                m[i][j] = truth_value_true;
            }
        }
    return m;
}

void MinimalityChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    adjacency_matrix_t m = matrix; // copy

    if (maximize) // flip all truth values
        m = flipMatrix(m);

    try
    {
        for (auto ordering : vertexOrderings)
        {
            auto matrixCopy = m;
            if (directed)
                checkMinimalityDir(matrixCopy, ordering, config);
            else
                checkMinimality(matrixCopy, ordering, config);
        }
    }
    catch (LimitReachedException e)
    {
        limit_reached++;
    }
    catch (minimalit_check_result_t e)
    {
        forbidden_graph_t signedEdges = e.clause;

        // in minimality check the oppposite signs were used, so have to flip it (for minimality)
        if (!maximize)
        {
            for (size_t i = 0; i < signedEdges.size(); i++)
                if (signedEdges[i].first == truth_value_true)
                    signedEdges[i].first = truth_value_false;
                else
                    signedEdges[i].first = truth_value_true;
        }

        if (symBreakClauses)
        {
            for (size_t i = 0; i < signedEdges.size(); i++)
                fprintf(symBreakClauses, "%s(%d,%d) ", signedEdges[i].first == truth_value_true ? "-" : "", signedEdges[i].second.first, signedEdges[i].second.second);
            // fprintf(symBreakClauses, "\n");

            fprintf(symBreakClauses, ";");
            // print permutation
            for (auto v : e.permutation)
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

        // TODO currently flipped edges but later replace it in the minimality check itself
        // flip sign of all the edges
        for (int m = 0; m < (int)matrices.size(); m++) // flip for all of them separately
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
