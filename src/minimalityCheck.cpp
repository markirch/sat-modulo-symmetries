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
    if (MAXIMIZE_SMS)
    {
        auto adj_copy = adjacency_matrix; // avoid side effects by copying
        // printf("Before:");
        // printAdjacencyMatrix(adj_copy);
        for (size_t i = 0; i < adjacency_matrix.size(); i++)
            for (size_t j = 0; j < adjacency_matrix.size(); j++)
            {
                if (adjacency_matrix[i][j] == truth_value_false)
                {
                    adj_copy[i][j] = truth_value_true;
                }
                else if (adjacency_matrix[i][j] == truth_value_true)
                {
                    adj_copy[i][j] = truth_value_false;
                }
            }
        // printf("Complement:");
        // printAdjacencyMatrix(adj_copy);
        int count = 0;
        isMinimal(vertex_ordering, config.initial_partition, 0, adj_copy, config, count);
    }
    else
    {
        int count = 0;
        isMinimal(vertex_ordering, config.initial_partition, 0, adjacency_matrix, config, count);
        // printf("Number of calls %d\n", count);
    }
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

        if (false) // check if for two in the partition vertices the neighborhood to the other vertices is the same
        {
            bool skip = false;
            int vertex = vertices[i];
            for (int j = row; j < i; j++)
            {
                int otherVertex = vertices[j];
                bool locallyEqual = true;
                for (int k = row; k < (int)vertices.size(); k++)
                {
                    if (k == i || k == j)
                        continue;
                    if (adjacency_matrix[vertex][vertices[k]] != adjacency_matrix[otherVertex][vertices[k]])
                    {
                        // printf("Witness: %d for %d %d\n", vertices[k], vertex, otherVertex);
                        locallyEqual = false;
                        break;
                    }
                }
                if (locallyEqual)
                {
                    skip = true;
                    break;
                }
            }
            if (skip)
            {
                continue;
            }
        }
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
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Depth: %d\n", row);
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("permutation: ");
    // for (auto v : vertices)
    //     printf("%d ", v);
    // printf("\n");
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Partition: ");
    // for (auto b : partition)
    //     if (b)
    //         printf("1");
    //     else
    //         printf("0");
    // printf("\n");
    vertex_t n = adjacency_matrix.size();
    vertex_t vertex = vertices[row];
    vertex_t col = row + 1; // current column

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
                    if ((vertex == row && vertices[col] == col) || (vertex == col && vertices[col] == row))
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
            // printf("row: %d; col: %d\n", row, i);
            // PRINT_CURRENT_LINE
            if (adjacency_matrix[row][i] == truth_value_true)
            {
                createClause(vertices, adjacency_matrix, config);
            }
            if (adjacency_matrix[row][i] == truth_value_false)
                return;

            partition[i] = 1;

            if (vertex == row)
            {
                if (vertices[i] == i) // usefull, if not permuted a lot yet
                    continue;

                int pos = getElementFromArray(&vertices[start], end - start, i); // eventually first search in unknown, but still need position for swapping.
                if (pos == -1)
                {
                    return;
                }
                std::swap(vertices[i], vertices[pos + start]);
            }
            else if (vertex == i)
            {
                if (vertices[i] == row)
                    continue;

                int pos = getElementFromArray(&vertices[start], end - start, row); // eventually first search in unknown, but still need position for swapping.
                if (pos == -1)
                    return;
                std::swap(vertices[i], vertices[pos + start]);
            }
            else
                return;
        }

        for (vertex_t i = start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i < end; i++) // check if remaining all are indeed adjacent
        {
            if (adjacency_matrix[row][i] != truth_value_true)
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
    /*
    printf("Add new Clause:\n");
    printAdjacencyMatrix(adjacency_matrix);
    for (int i = 0; i < adjacency_matrix.size(); i++)
        printf("%d ", vertices[i]);
    printf("\n"); */

    int n = adjacency_matrix.size();
    std::vector<signed_edge_t> edges;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if ((i == vertices[i] && j == vertices[j]) || (i == vertices[j] && j == vertices[i]))
                continue;

            int isAdjacentNormal = adjacency_matrix[i][j];
            int isAdjacentPerm = adjacency_matrix[vertices[i]][vertices[j]];

            if ((isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_unknown))
            {
                // printf("faulty %d %d:\n", i, j);
                // printAdjacencyMatrix(adjacency_matrix);

                // printf("Permutation: ");
                // for (int i = 0; i < adjacency_matrix.size(); i++)
                //     printf("%d ", vertices[i]);
                // printf("\n");

                // adjacency_matrix_t permAdj = adjacency_matrix;
                // for (int i = 0; i < adjacency_matrix.size(); i++)
                //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
                //     {
                //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
                //     }
                // printf("Permuted matrix:\n");
                // printAdjacencyMatrix(permAdj);

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
                //     for (int i = 0; i < (int)adjacency_matrix.size(); i++)
                //         fprintf(symBreakClausesFile, "%d ", vertices[i]);
                //     fprintf(symBreakClausesFile, "; ");
                // }
                if (MAXIMIZE_SMS)
                {
                    // flip sign of all the edges
                    for (size_t i = 0; i < edges.size(); i++)
                        if (edges[i].first == truth_value_true)
                            edges[i].first = truth_value_false;
                        else
                            edges[i].first = truth_value_true;
                }
                minimalit_check_result_t res;
                res.permutation = vertices;
                res.clause = edges;
                throw res;
            }

            if (isAdjacentNormal == truth_value_true)
                edges.push_back({truth_value_false, {i, j}});
            else
                edges.push_back({truth_value_true, {vertices[i], vertices[j]}});
        }
    }
    // printf("Add new Clause:\n");
    // printAdjacencyMatrix(adjacency_matrix);
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     printf("%d ", vertices[i]);
    // printf("\n");

    // adjacency_matrix_t permAdj = adjacency_matrix;
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
    //     {
    //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
    //     }
    // printAdjacencyMatrix(permAdj);

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

// ------------------------------------ version with multiple graphs over the same vertex set. Minimal by the following ordering: first row of first graph, first row of second graph first row of third graph,..., second row of first graph,...

void createClauseMultiple(vertex_ordering_t &vertices, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t);
void isMinimalVertexMultiple(vertex_ordering_t vertices, partition_t partition, int row, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t config, int &count);
void isMinimalMultiple(vertex_ordering_t vertices, partition_t partition, int row, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t config, int &count);

void checkMinimalityMultiple(vector<adjacency_matrix_t> &adjacency_matrices, vertex_ordering_t vertex_ordering, minimalit_check_config_t config)
{

    int count = 0;
    isMinimalMultiple(vertex_ordering, config.initial_partition, 0, adjacency_matrices, config, count);
    // printf("Number of calls %d\n", count);
}

void isMinimalMultiple(vertex_ordering_t vertices, partition_t partition, int row, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t config, int &count)
{
    // PRINT_CURRENT_LINE
    count++;
    if (config.cutoff != 0 && count > config.cutoff)
    {

        reachedCutoff++;
        throw LimitReachedException();
    }

    int n = vertices.size();
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

        isMinimalVertexMultiple(newVertices, newPartition, row, adjacency_matrices, config, count);
    }
}

// Adapt the partition according the selected vertex for the row
void isMinimalVertexMultiple(vertex_ordering_t vertices, partition_t partition, int row, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t config, int &count)
{
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Depth: %d\n", row);
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("permutation: ");
    // for (auto v : vertices)
    //     printf("%d ", v);
    // printf("\n");
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Partition: ");
    // for (auto b : partition)
    //     if (b)
    //         printf("1");
    //     else
    //         printf("0");
    // printf("\n");
    for (adjacency_matrix_t adjacency_matrix : adjacency_matrices)
    {
        vertex_t n = vertices.size();
        vertex_t vertex = vertices[row];
        vertex_t col = row + 1; // current column

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
                        if ((vertex == row && vertices[col] == col) || (vertex == col && vertices[col] == row))
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
                    createClauseMultiple(vertices, adjacency_matrices, config);
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
                    createClauseMultiple(vertices, adjacency_matrices, config);
            }

            if (start + (vertex_t)notAdjacentList.size() < end && adjacency_matrix[row][start + notAdjacentList.size()] == truth_value_false)
                return;

            // match undefined
            for (vertex_t i = start + notAdjacentList.size(); i < start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i++)
            {
                // printf("row: %d; col: %d\n", row, i);
                // PRINT_CURRENT_LINE
                if (adjacency_matrix[row][i] == truth_value_true)
                {
                    createClauseMultiple(vertices, adjacency_matrices, config);
                }
                if (adjacency_matrix[row][i] == truth_value_false)
                    return;

                partition[i] = 1;

                if (vertex == row)
                {
                    if (vertices[i] == i) // usefull, if not permuted a lot yet
                        continue;

                    int pos = getElementFromArray(&vertices[start], end - start, i); // eventually first search in unknown, but still need position for swapping.
                    if (pos == -1)
                    {
                        return;
                    }
                    std::swap(vertices[i], vertices[pos + start]);
                }
                else if (vertex == i)
                {
                    if (vertices[i] == row)
                        continue;

                    int pos = getElementFromArray(&vertices[start], end - start, row); // eventually first search in unknown, but still need position for swapping.
                    if (pos == -1)
                        return;
                    std::swap(vertices[i], vertices[pos + start]);
                }
                else
                    return;
            }

            for (vertex_t i = start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i < end; i++) // check if remaining all are indeed adjacent
            {
                if (adjacency_matrix[row][i] != truth_value_true)
                    return;
            }

            partition[start + notAdjacentList.size() + unknownList.size()] = 1; // start of partition with one list
            col = end;
        }
    }

    return isMinimalMultiple(vertices, partition, row + 1, adjacency_matrices, config, count);
}

void createClauseMultiple(vertex_ordering_t &vertices, vector<adjacency_matrix_t> &adjacency_matrices, minimalit_check_config_t)
{
    /*
    printf("Add new Clause:\n");
    printAdjacencyMatrix(adjacency_matrix);
    for (int i = 0; i < adjacency_matrix.size(); i++)
        printf("%d ", vertices[i]);
    printf("\n"); */

    int n = (int)vertices.size();
    int m = (int)adjacency_matrices.size();
    vector<vector<signed_edge_t>> edgesMulti(m); // for each graph generate a list for the edges

    for (int i = 0; i < n; i++)
    {
        for (int a = 0; a < m; a++)
        {
            adjacency_matrix_t adjacency_matrix = adjacency_matrices[a];
            for (int j = i + 1; j < n; j++)
            {
                if ((i == vertices[i] && j == vertices[j]) || (i == vertices[j] && j == vertices[i]))
                    continue;

                int isAdjacentNormal = adjacency_matrix[i][j];
                int isAdjacentPerm = adjacency_matrix[vertices[i]][vertices[j]];

                if ((isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_true) ||
                    (isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_unknown) ||
                    (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_true) ||
                    (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_unknown))
                {
                    // printf("faulty %d %d:\n", i, j);
                    // printAdjacencyMatrix(adjacency_matrix);

                    // printf("Permutation: ");
                    // for (int i = 0; i < adjacency_matrix.size(); i++)
                    //     printf("%d ", vertices[i]);
                    // printf("\n");

                    // adjacency_matrix_t permAdj = adjacency_matrix;
                    // for (int i = 0; i < adjacency_matrix.size(); i++)
                    //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
                    //     {
                    //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
                    //     }
                    // printf("Permuted matrix:\n");
                    // printAdjacencyMatrix(permAdj);

                    EXIT_UNWANTED_STATE
                }

                if ((isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_false) ||
                    (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_unknown) ||
                    (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_false))
                {
                    edgesMulti[a].push_back(std::make_pair(truth_value_false, std::make_pair(i, j)));
                    edgesMulti[a].push_back(std::make_pair(truth_value_true, std::make_pair(vertices[i], vertices[j])));

                    // if (symBreakClausesFile)
                    // {
                    //     for (int i = 0; i < (int)adjacency_matrix.size(); i++)
                    //         fprintf(symBreakClausesFile, "%d ", vertices[i]);
                    //     fprintf(symBreakClausesFile, "; ");
                    // }
                    if (MAXIMIZE_SMS)
                    {
                        // flip sign of all the edges
                        for (int a = 0; a < m; a++)
                            for (size_t i = 0; i < edgesMulti[a].size(); i++)
                                if (edgesMulti[a][i].first == truth_value_true)
                                    edgesMulti[a][i].first = truth_value_false;
                                else
                                    edgesMulti[a][i].first = truth_value_true;
                    }
                    minimalit_check_result_multi_t res;
                    res.permutation = vertices;
                    res.clause = edgesMulti;
                    throw res;
                }

                if (isAdjacentNormal == truth_value_true)
                    edgesMulti[a].push_back({truth_value_false, {i, j}});
                else
                    edgesMulti[a].push_back({truth_value_true, {vertices[i], vertices[j]}});
            }
        }
    }
    // printf("Add new Clause:\n");
    // printAdjacencyMatrix(adjacency_matrix);
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     printf("%d ", vertices[i]);
    // printf("\n");

    // adjacency_matrix_t permAdj = adjacency_matrix;
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
    //     {
    //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
    //     }
    // printAdjacencyMatrix(permAdj);

    EXIT_UNWANTED_STATE
}

// -----------------------------complementVersion-------------------------------------
void createClauseComplement(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t);
void isMinimalComplement(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count);
void isMinimalVertexComplement(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count);
void createClauseComplement(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t);

static inline truth_value_t invertState(truth_value_t state)
{
    if (state == truth_value_true)
        return truth_value_false;
    if (state == truth_value_false)
        return truth_value_true;
    return truth_value_unknown;
}

void checkMinimalityComplement(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config)
{
    int count = 0;
    isMinimalComplement(vertex_ordering, config.initial_partition, 0, adjacency_matrix, config, count);
    // printf("Number of calls %d\n", count);
}

void isMinimalComplement(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count)
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

        isMinimalVertexComplement(newVertices, newPartition, row, adjacency_matrix, config, count);
    }
}

// Adapt the partition according the selected vertex for the row
void isMinimalVertexComplement(vertex_ordering_t vertices, partition_t partition, int row, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t config, int &count)
{
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Depth: %d\n", row);
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("permutation: ");
    // for (auto v : vertices)
    //     printf("%d ", v);
    // printf("\n");
    // for (int i = 0; i < row; i++)
    //     printf("   ");
    // printf("Partition: ");
    // for (auto b : partition)
    //     if (b)
    //         printf("1");
    //     else
    //         printf("0");
    // printf("\n");
    vertex_t n = adjacency_matrix.size();
    vertex_t vertex = vertices[row];
    vertex_t col = row + 1; // current column

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
            statePermutation = invertState(statePermutation);

            if (state == statePermutation)
            {
                if (state != truth_value_unknown)
                {
                    col = end;
                    continue;
                }
                else
                {
                    if ((vertex == row && vertices[col] == col) || (vertex == col && vertices[col] == row))
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
                createClauseComplement(vertices, adjacency_matrix, config);
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
            int state = invertState(adjacency_matrix[vertex][vertices[i]]);

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
                createClauseComplement(vertices, adjacency_matrix, config);
        }

        if (start + (vertex_t)notAdjacentList.size() < end && adjacency_matrix[row][start + notAdjacentList.size()] == truth_value_false)
            return;

        // match undefined
        for (vertex_t i = start + notAdjacentList.size(); i < start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i++)
        {
            // printf("row: %d; col: %d\n", row, i);
            // PRINT_CURRENT_LINE
            if (adjacency_matrix[row][i] == truth_value_true)
            {
                createClauseComplement(vertices, adjacency_matrix, config);
            }
            if (adjacency_matrix[row][i] == truth_value_false)
                return;

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
                createClauseComplement(vertices, adjacency_matrix, config); // original must be 0 because then other is 1
            }
            else if (vertex == i)
            {
                if (vertices[i] == row)
                    continue;

                int pos = getElementFromArray(&vertices[start], end - start, row); // eventually first search in truth_value_unknown, but still need position for swapping.
                if (pos == -1)
                    return;
                std::swap(vertices[i], vertices[pos + start]);
                createClauseComplement(vertices, adjacency_matrix, config); // original must be 0 because then other is 1
            }
            else
                return;
        }

        for (vertex_t i = start + (vertex_t)(notAdjacentList.size() + unknownList.size()); i < end; i++) // check if remaining all are indeed adjacent
        {
            if (adjacency_matrix[row][i] != truth_value_true)
                return;
        }

        partition[start + notAdjacentList.size() + unknownList.size()] = 1; // start of partition with one list
        col = end;
    }

    return isMinimalComplement(vertices, partition, row + 1, adjacency_matrix, config, count);
}

void createClauseComplement(vertex_ordering_t &vertices, adjacency_matrix_t &adjacency_matrix, minimalit_check_config_t)
{

    int n = adjacency_matrix.size();
    std::vector<signed_edge_t> edges;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            // equality of edges is not considered anymore

            int isAdjacentNormal = adjacency_matrix[i][j];
            int isAdjacentPerm = invertState(adjacency_matrix[vertices[i]][vertices[j]]);

            if ((i == vertices[i] && j == vertices[j]) || (i == vertices[j] && j == vertices[i]))
            {
                if (isAdjacentNormal != truth_value_false)
                {
                    edges.push_back({truth_value_false, {i, j}});
                    minimalit_check_result_t res;
                    res.permutation = vertices;
                    res.clause = edges;
                    throw res;
                }
                return;
            }

            if ((isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_false && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_true) ||
                (isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_unknown))
            {
                // printf("faulty %d %d:\n", i, j);
                // printAdjacencyMatrix(adjacency_matrix);

                // printf("Permutation: ");
                // for (int i = 0; i < adjacency_matrix.size(); i++)
                //     printf("%d ", vertices[i]);
                // printf("\n");

                // adjacency_matrix_t permAdj = adjacency_matrix;
                // for (int i = 0; i < adjacency_matrix.size(); i++)
                //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
                //     {
                //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
                //     }
                // printf("Permuted matrix:\n");
                // printAdjacencyMatrix(permAdj);

                EXIT_UNWANTED_STATE
            }

            if ((isAdjacentNormal == truth_value_unknown && isAdjacentPerm == truth_value_false) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_unknown) ||
                (isAdjacentNormal == truth_value_true && isAdjacentPerm == truth_value_false))
            {
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(i, j)));
                edges.push_back(std::make_pair(truth_value_false, std::make_pair(vertices[i], vertices[j])));

                // if (symBreakClausesFile)
                // {
                //     for (int i = 0; i < (int)adjacency_matrix.size(); i++)
                //         fprintf(symBreakClausesFile, "%d ", vertices[i]);
                //     fprintf(symBreakClausesFile, "; ");
                // }
                minimalit_check_result_t res;
                res.permutation = vertices;
                res.clause = edges;
                throw res;
            }

            if (isAdjacentNormal == truth_value_true)
                edges.push_back({truth_value_false, {i, j}});
            else
                edges.push_back({truth_value_false, {vertices[i], vertices[j]}});
        }
    }
    // printf("Add new Clause:\n");
    // printAdjacencyMatrix(adjacency_matrix);
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     printf("%d ", vertices[i]);
    // printf("\n");

    // adjacency_matrix_t permAdj = adjacency_matrix;
    // for (int i = 0; i < adjacency_matrix.size(); i++)
    //     for (int j = i + 1; j < adjacency_matrix.size(); j++)
    //     {
    //         permAdj[j][i] = permAdj[i][j] = adjacency_matrix[vertices[i]][vertices[j]];
    //     }
    // printAdjacencyMatrix(permAdj);

    EXIT_UNWANTED_STATE
}
