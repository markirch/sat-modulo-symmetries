#include "minimalityCheck.hpp"

// Set to 1 for colex maximality, 0 for colex minimality.
#define COLEX_MAXIMIZE 0

static inline bool permutation_is_strictly_better(truth_value_t state, truth_value_t statePerm)
{
#if COLEX_MAXIMIZE
    return (state == truth_value_unknown && statePerm == truth_value_true) ||
           (state == truth_value_false && statePerm == truth_value_unknown) ||
           (state == truth_value_false && statePerm == truth_value_true);
#else
    return (state == truth_value_true && statePerm == truth_value_unknown) ||
           (state == truth_value_true && statePerm == truth_value_false) ||
           (state == truth_value_unknown && statePerm == truth_value_false);
#endif
}

static inline bool permutation_cannot_improve(truth_value_t state, truth_value_t statePerm)
{
#if COLEX_MAXIMIZE
    return (state == truth_value_true && statePerm == truth_value_unknown) ||
           (state == truth_value_true && statePerm == truth_value_false) ||
           (state == truth_value_unknown && statePerm == truth_value_false) ||
           (state == truth_value_unknown && statePerm == truth_value_unknown);
#else
    return (state == truth_value_false && statePerm == truth_value_unknown) ||
           (state == truth_value_false && statePerm == truth_value_true) ||
           (state == truth_value_unknown && statePerm == truth_value_true) ||
           (state == truth_value_unknown && statePerm == truth_value_unknown);
#endif
}

// internal functions
void isMinimalColex(vertex_ordering_t vertices, int col, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count);
void createClauseColex(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config);

void checkMinimalityColex(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config)
{
    PRINT_CURRENT_LINE

    int count = 0;
    isMinimalColex(vertex_ordering, 0, adjacency_matrix, config, count);
    // printf("Number of calls %d\n", count);
}

void isMinimalColex(vertex_ordering_t vertices, int col, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count)
{
    count++;
    if (config.cutoff != 0 && count > config.cutoff)
    {
        throw LimitReachedException();
    }

    int n = (int)adjacency_matrix.size();
    if (col == n)
        return;

    for (int i = col; i < n; i++)
    {
        std::swap(vertices[col], vertices[i]);

        // Compare only the current column (colex order).
        bool shouldRecurse = true;
        for (int j = 0; j < col; j++)
        {
            truth_value_t state = adjacency_matrix[j][col];
            truth_value_t statePerm = adjacency_matrix[vertices[j]][vertices[col]];

            if ((j == vertices[j] && col == vertices[col]) ||
                (j == vertices[col] && col == vertices[j]))
                continue;

            if (permutation_is_strictly_better(state, statePerm))
            {
                createClauseColex(vertices, adjacency_matrix, config);
            }

            if (permutation_cannot_improve(state, statePerm))
            {
                shouldRecurse = false;
                break;
            }

            // if ((state == truth_value_unknown) || (statePerm == truth_value_unknown))
            // {
            //     shouldRecurse = false;
            //     break;
            // }

            // all other cases are inconclusive under partial assignments
        }

        if (shouldRecurse)
            isMinimalColex(vertices, col + 1, adjacency_matrix, config, count);

        std::swap(vertices[col], vertices[i]);
    }
}

void createClauseColex(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t)
{
    int n = (int)adjacency_matrix.size();
    std::vector<signed_edge_t> edges;

    for (int j = 1; j < n; j++) // column
    {
        for (int i = 0; i < j; i++) // row
        {
            if ((i == vertices[i] && j == vertices[j]) || (i == vertices[j] && j == vertices[i]))
                continue;

            int isAdjacentNormal = adjacency_matrix[i][j];
            int isAdjacentPerm = adjacency_matrix[vertices[i]][vertices[j]];

            if (permutation_cannot_improve((truth_value_t)isAdjacentNormal, (truth_value_t)isAdjacentPerm))
            {
                // no clause can be created by the given permutation
                EXIT_UNWANTED_STATE
            }

            if (permutation_is_strictly_better((truth_value_t)isAdjacentNormal, (truth_value_t)isAdjacentPerm))
            {
#if COLEX_MAXIMIZE
                edges.push_back(std::make_pair(truth_value_true, std::make_pair(i, j)));
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(vertices[i], vertices[j])));
#else
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(i, j)));
                edges.push_back(std::make_pair(truth_value_true, std::make_pair(vertices[i], vertices[j])));
#endif

                minimalit_check_result_t res;
                res.permutation = vertices;
                res.clause = edges;
                throw res;
            }

            // Keep all previous comparisons equal so the first strict improvement remains decisive.
#if COLEX_MAXIMIZE
            if (isAdjacentNormal == truth_value_false)
                edges.push_back({truth_value_true, {i, j}});
            else
                edges.push_back({truth_value_false, {vertices[i], vertices[j]}});
#else
            if (isAdjacentNormal == truth_value_true)
                edges.push_back({truth_value_false, {i, j}});
            else
                edges.push_back({truth_value_true, {vertices[i], vertices[j]}});
#endif
        }
    }

    EXIT_UNWANTED_STATE
}
