#include "minimalityCheck.hpp"


int reachedCutoff = 0;

// internal functions
void isMinimal(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count);
void isMinimalVertex(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count);
void createClause(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config);

// TODO replace with vector function
int getElementFromArray(int *array, int n, int elem);
// void checkIfPartition(vertex_ordering_t &vertices);

void checkMinimality(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config)
{
    // printf("min check\n");
    int count = 0;
    if (MAXIMIZE_SMS)
    {
        adjacency_matrix_t adjacency_matrix_copy = adjacency_matrix;
        // flip adjacent to non adjacent and vice versa
        for (int i = 0; i < (int) adjacency_matrix.size(); i++)
        {
            for (int j = 0; j < (int) adjacency_matrix.size(); j++)
            {
                if (adjacency_matrix[i][j] == truth_value_true)
                    adjacency_matrix_copy[i][j] = truth_value_false;
                else if (adjacency_matrix[i][j] == truth_value_false)
                    adjacency_matrix_copy[i][j] = truth_value_true;
            }
        }
        isMinimal(vertex_ordering, config.initial_partition, 0, adjacency_matrix_copy, config, count);
    } else {
        isMinimal(vertex_ordering, config.initial_partition, 0, adjacency_matrix, config, count);
    }
    // printf("Number of calls %d\n", count);
}

void isMinimal(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count)
{
    // PRINT_CURRENT_LINE
    count++;
    if (config.cutoff != 0 && count > config.cutoff)
    {

        reachedCutoff++;
        throw LimitReachedException();
    }

    int n = adjacency_matrix.size();
    if (row == n)
        return;

    for (int i = row; i < n && (partition[i] == 0 || i == row); i++) // iter over vertices from partition.
    {

        // select vertices[i] to be on position row

        partition_t newPartition = partition; // copy partition (works because vectors)
        vertex_ordering_t newVertices = vertices;

        newPartition[row + 1] = 1;
        newVertices[row] = vertices[i];
        newVertices[i] = vertices[row];

        isMinimalVertex(newVertices, newPartition, row, adjacency_matrix, config, count);
    }
}

// Adapt the partition according the selected vertex for the row
void isMinimalVertex(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count)
{
    vertex_t n = adjacency_matrix.size();
    vertex_t vertex = vertices[row];
    vertex_t col = row + 1; // current column

    /*for (int i =0 ; i < row; i++)
        printf("\t");
    printf("%d\n", vertex); */

    while (col != n) // iter over all partitions
    {

        vertex_t start = col;
        vertex_t end;
        for (end = col + 1; end < n && partition[end] == 0; end++)
            ; // get the single partition

        if (start + 1 == end) // single vertex partition as seperate case (minor improvement)
        {

            truth_value_t state, statePermutation;
            state = adjacency_matrix[row][start];
            statePermutation = adjacency_matrix[vertex][vertices[start]];

            if (state == statePermutation)
            {
                if (state != truth_value_unknown)
                {
                    col = end;
                    continue;
                }
                else
                {
                    if (vertex == row && vertices[col] == col)
                    {
                        col = end;
                        continue;
                    }
                }
                return;
            }

            if ((state == truth_value_unknown && statePermutation == truth_value_false) ||
                (state == truth_value_true && statePermutation == truth_value_unknown) ||
                (state == truth_value_true && statePermutation == truth_value_false))
            {
                createClause(vertices, adjacency_matrix, config);
                col = end;
                continue;
            }

            return;
        }

        std::vector<vertex_t> adjacentList(0);
        std::vector<vertex_t> unknownList(0);
        std::vector<vertex_t> notAdjacentList(0);
        for (int i = start; i < end; i++)
        {
            int state = adjacency_matrix[vertex][vertices[i]];

            if (state == truth_value_unknown)
                unknownList.push_back(vertices[i]);
            else if (state == truth_value_true)
                adjacentList.push_back(vertices[i]);
            else if (state == truth_value_false)
                notAdjacentList.push_back(vertices[i]);
        }

        if (start + (vertex_t)notAdjacentList.size() < n) // create new partition containing all the non adjacent once.
            partition[start + notAdjacentList.size()] = 1;

        // printf("%d,%d,%d,%d,%d,%d\n", start, nNotAdjacent, start + nNotAdjacent, nUnknown, start + nNotAdjacent + nUnknown, nAdjacent);
        std::copy(notAdjacentList.begin(), notAdjacentList.end(), &vertices[start]);
        std::copy(unknownList.begin(), unknownList.end(), &vertices[start + notAdjacentList.size()]);
        std::copy(adjacentList.begin(), adjacentList.end(), &vertices[start + notAdjacentList.size() + unknownList.size()]);

        // iter over zeros and check if at least as many zeros
        for (vertex_t i = start; i < start + (vertex_t)notAdjacentList.size(); i++)
        {
            if (adjacency_matrix[row][i] != truth_value_false)
                createClause(vertices, adjacency_matrix, config);
        }

        if (start + (vertex_t)notAdjacentList.size() < end && adjacency_matrix[row][start + notAdjacentList.size()] == truth_value_false)
            return;

        // match undefined
        for (vertex_t i = start + notAdjacentList.size(); i < start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i++)
        {
            if (adjacency_matrix[row][i] == truth_value_true)
                createClause(vertices, adjacency_matrix, config);

            partition[i] = 1;

            if (vertex == row)
            {
                if (vertices[i] == i) // usefull, if not permuted a lot yet
                    continue;

                int pos = getElementFromArray(&vertices[start], end - start, i); // eventually first search in truth_value_unknown, but still need position for swapping.
                if (pos == -1)
                {
                    return;
                }
                std::swap(vertices[i], vertices[pos + start]);
                continue;
            }
            return;
        }

        if ((vertex_t)(start + notAdjacentList.size() + unknownList.size()) < end && adjacency_matrix[row][start + notAdjacentList.size() + unknownList.size()] == truth_value_unknown)
            return;

        for (vertex_t i = start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i < end; i++)
        {
            if (adjacency_matrix[row][i] != truth_value_true)
                return;
        }

        partition[start + notAdjacentList.size() + unknownList.size()] = 1; // start of partition with one list
        col = end;
    }

    // PRINT_CURRENT_LINE

    // the same for incomming edges TODO refactor
    col = row + 1;
    while (col != n) // iter over all partitions
    {

        vertex_t start = col;
        vertex_t end;
        for (end = col + 1; end < n && partition[end] == 0; end++)
            ; // get the single partition

        if (start + 1 == end) // single vertex partition as seperate case (minor improvement)
        {

            truth_value_t state, statePermutation;
            state = adjacency_matrix[start][row];
            statePermutation = adjacency_matrix[vertices[start]][vertex];

            if (state == statePermutation)
            {
                if (state != truth_value_unknown)
                {
                    col = end;
                    continue;
                }
                else
                {
                    if (vertex == row && vertices[col] == col)
                    {
                        col = end;
                        continue;
                    }
                }
                return;
            }

            if ((state == truth_value_unknown && statePermutation == truth_value_false) ||
                (state == truth_value_true && statePermutation == truth_value_unknown) ||
                (state == truth_value_true && statePermutation == truth_value_false))
            {
                createClause(vertices, adjacency_matrix, config);
                col = end;
                continue;
            }

            return;
        }

        std::vector<vertex_t> adjacentList(0);
        std::vector<vertex_t> unknownList(0);
        std::vector<vertex_t> notAdjacentList(0);
        for (int i = start; i < end; i++)
        {
            int state = adjacency_matrix[vertices[i]][vertex];

            if (state == truth_value_unknown)
                unknownList.push_back(vertices[i]);
            else if (state == truth_value_true)
                adjacentList.push_back(vertices[i]);
            else if (state == truth_value_false)
                notAdjacentList.push_back(vertices[i]);
        }

        if (start + (vertex_t)notAdjacentList.size() < n) // create new partition containing all the non adjacent once.
            partition[start + notAdjacentList.size()] = 1;

        // printf("%d,%d,%d,%d,%d,%d\n", start, nNotAdjacent, start + nNotAdjacent, nUnknown, start + nNotAdjacent + nUnknown, nAdjacent);
        std::copy(notAdjacentList.begin(), notAdjacentList.end(), &vertices[start]);
        std::copy(unknownList.begin(), unknownList.end(), &vertices[start + notAdjacentList.size()]);
        std::copy(adjacentList.begin(), adjacentList.end(), &vertices[start + notAdjacentList.size() + unknownList.size()]);

        // iter over zeros and check if at least as many zeros
        for (vertex_t i = start; i < start + (vertex_t)notAdjacentList.size(); i++)
        {
            if (adjacency_matrix[i][row] != truth_value_false)
                createClause(vertices, adjacency_matrix, config);
        }

        if (start + (vertex_t)notAdjacentList.size() < end && adjacency_matrix[start + notAdjacentList.size()][row] == truth_value_false)
            return;

        // match undefined
        for (vertex_t i = start + notAdjacentList.size(); i < start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i++)
        {
            if (adjacency_matrix[i][row] == truth_value_true)
                createClause(vertices, adjacency_matrix, config);

            partition[i] = 1;

            if (vertex == row)
            {
                if (vertices[i] == i) // usefull, if not permuted a lot yet
                    continue;

                int pos = getElementFromArray(&vertices[start], end - start, i); // eventually first search in truth_value_unknown, but still need position for swapping.
                if (pos == -1)
                {
                    return;
                }
                std::swap(vertices[i], vertices[pos + start]);
                continue;
            }
            return;
        }

        if ((vertex_t)(start + notAdjacentList.size() + unknownList.size()) < end && adjacency_matrix[start + notAdjacentList.size() + unknownList.size()][row] == truth_value_unknown)
            return;

        for (vertex_t i = start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i < end; i++)
        {
            if (adjacency_matrix[i][row] != truth_value_true)
                return;
        }

        partition[start + notAdjacentList.size() + unknownList.size()] = 1; // start of partition with one list
        col = end;
    }

    return isMinimal(vertices, partition, row + 1, adjacency_matrix, config, count);
}

// TODO replace with vector function
int getElementFromArray(int *array, int n, int elem)
{
    for (int i = 0; i < n; i++)
    {
        if (array[i] == elem)
            return i;
    }
    return -1;
}

void createClause(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t)
{
    // printf("create clause\n");
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //    printf("%d ", vertices[i]);
    // printf("\n");

    int n = adjacency_matrix.size();
    std::vector<signed_edge_t> edges;

    for (int i = 0; i < n; i++)
    {
        // outgoing edges
        for (int j = i + 1; j < n; j++)
        {
            if (i == vertices[i] && j == vertices[j]) // skip if the same **directed** edge
                continue;

            int isAdjacentNormal = adjacency_matrix[i][j];
            int isAdjacentPerm = adjacency_matrix[vertices[i]][vertices[j]];

            if ((isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_true))
            {
                EXIT_UNWANTED_STATE
            }

            if ((isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_false) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_false))
            {
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(i, j)));
                edges.push_back(std::make_pair(truth_value_true, std::make_pair(vertices[i], vertices[j])));

                // if (symBreakClausesFile)
                // {
                //     for (int i = 0; i < adjacency_matrix.size(); i++)
                //         fprintf(symBreakClausesFile, "%d ", vertices[i]);
                //     fprintf(symBreakClausesFile, "; ");
                // }
                minimalit_check_result_t res;
                res.permutation = vertices;
                res.clause = edges;
                if (MAXIMIZE_SMS)
                {
                    // reverse signes of edges
                    for (int i = 0; i < (int)edges.size(); i++)
                    {
                        edges[i].first = (edges[i].first == truth_value_true) ? truth_value_false : truth_value_true;
                    }
                }
                throw res;
            }

            if (isAdjacentNormal == truth_value_true)
                edges.push_back({truth_value_false, {i, j}});
            else
                edges.push_back({truth_value_true, {vertices[i], vertices[j]}});
        }

        // incomming edges TODO refactor because almost the same as above (just swapped i and j for the edges)
        for (int j = i + 1; j < n; j++)
        {
            if (i == vertices[i] && j == vertices[j])
                continue;

            int isAdjacentNormal = adjacency_matrix[j][i];
            int isAdjacentPerm = adjacency_matrix[vertices[j]][vertices[i]];

            if ((isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_true))
            {
                return;
            }

            if ((isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_false) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_false))
            {
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(j, i)));
                edges.push_back(std::make_pair(truth_value_true, std::make_pair(vertices[j], vertices[i])));

                // if (symBreakClausesFile)
                // {
                //     for (int i = 0; i < adjacency_matrix.size(); i++)
                //         fprintf(symBreakClausesFile, "%d ", vertices[i]);
                //     fprintf(symBreakClausesFile, "; ");
                // }
                minimalit_check_result_t res;
                res.permutation = vertices;
                res.clause = edges;
                throw res;
            }

            if (isAdjacentNormal == truth_value_true)
                edges.push_back({truth_value_false, {j, i}});
            else
                edges.push_back({truth_value_true, {vertices[j], vertices[i]}});
        }
    }
    EXIT_UNWANTED_STATE
}

// ------------------Functions for debugging ----------------------------------

// check if all vertices are contained
// void checkIfPartition(std::vector<vertex_t> vertices)
// {
//     std::vector<bool> reached(vertices.size(), false);
//     for (vertex_t i = 0; i < (vertex_t)vertices.size(); i++)
//         reached[vertices[i]] = true; // mark each found vertex to true

//     if (std::any_of(reached.begin(), reached.end(), [](auto x)
//                     { return !x; })) // check if no doubles
//         EXIT_UNWANTED_STATE
// }

void checkMinimalityMultiple(vector<adjacency_matrix_t> &, vertex_ordering_t, minimalit_check_config_t)
{
    EXIT_UNWANTED_STATE // not implemented yet
}