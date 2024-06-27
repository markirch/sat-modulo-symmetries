#ifndef MIN_CHECK_H
#define MIN_CHECK_H

#include "useful.h"
#include "graphChecker.hpp"

#define MAXIMIZE_SMS 0

class LimitReachedException
{
};

typedef struct
{
    vector<int> permutation;
    vector<signed_edge_t> clause; // at least one of this edges/non-edges must be present
} minimalit_check_result_t;

typedef struct
{
    vector<int> permutation;
    vector<vector<signed_edge_t>> clause; // at least one of this edges/non-edges must be present in one of the graphs
} minimalit_check_result_multi_t;

typedef std::vector<vertex_t> vertex_ordering_t;
typedef std::vector<bool> partition_t; // if true, then new partition starts

typedef struct
{
    partition_t initial_partition; // true if new partition starts
    int cutoff;
} minimalit_check_config_t; // used for different configurations of the minimality check

/**
 * Check if graph is minimal and add a clause if this is not the case.
 * Throws std::vector<signed_edge_t> if current partially defined graph is not minimal
 * Throws LimitReachedException if cutoff limit is reached
 */
void checkMinimality(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityComplement(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityMultiple(vector<adjacency_matrix_t> &adjacency_matrices, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);

// basic minimality check with potentially different initial vertexOrderings
class MinimalityChecker : public PartiallyDefinedGraphChecker
{
    vector<vertex_ordering_t> vertexOrderings;
    minimalit_check_config_t config;
    FILE *symBreakClauses;
public:
    MinimalityChecker(int frequency, partition_t initial_partition, vector<vertex_ordering_t> vertexOrderings, int cutoff, FILE *symBreakClauses)
    {
        this->name = "MinimalityChecker";
        this->frequency = frequency;
        this->config.initial_partition = initial_partition;
        this->vertexOrderings = vertexOrderings;
        this->config.cutoff = cutoff;
        this->symBreakClauses = symBreakClauses;
    }

    void checkProperty(const adjacency_matrix_t &matrix);
};

// basic minimality check with potentially different initial vertexOrderings
class MultipleMinimalityChecker : public PartiallyDefinedMultiGraphChecker
{
    vector<vertex_ordering_t> vertexOrderings;
    minimalit_check_config_t config;
    FILE *symBreakClauses;
public:
    MultipleMinimalityChecker(int frequency, partition_t initial_partition, vector<vertex_ordering_t> vertexOrderings, int cutoff, FILE *symBreakClauses)
    {
        this->name = "MinimalityCheckerMulti";
        this->frequency = frequency;
        this->config.initial_partition = initial_partition;
        this->vertexOrderings = vertexOrderings;
        this->config.cutoff = cutoff;
        this->symBreakClauses = symBreakClauses;
    }

    void checkProperty(const vector<adjacency_matrix_t> &matrices);
};

// basic minimality check with potentially different initial vertexOrderings
class MaximalityChecker : public PartiallyDefinedGraphChecker
{
    vector<vertex_ordering_t> vertexOrderings;
    minimalit_check_config_t config;

public:
    MaximalityChecker(int frequency, partition_t initial_partition, vector<vertex_ordering_t> vertexOrderings, int cutoff)
    {
        this->name = "MaximalityChecker";
        this->frequency = frequency;
        this->config.initial_partition = initial_partition;
        this->vertexOrderings = vertexOrderings;
        this->config.cutoff = cutoff;
    }

    void checkProperty(const adjacency_matrix_t &matrix);
};

// check minimality based on preordering of vertices given by the truth values of certain variables
class MinimalityCheckerWithStaticPartition : public ComplexPartiallyDefinedGraphChecker
{
    int cutoff;
    edge_vars_t edges;
    edge_vars_t staticParitionVars; // variables describing the static partition.

public:
    MinimalityCheckerWithStaticPartition(int frequency, int cutoff, edge_vars_t edges, edge_vars_t staticParitionVars)
    {
        this->frequency = frequency;
        this->cutoff = cutoff;
        this->edges = edges;
        this->staticParitionVars = staticParitionVars;
    }

    void checkProperty(const adjacency_matrix_t &matrix, const vector<truth_value_t> &currentAssignemnt); // { EXIT_UNWANTED_STATE }
};

/* TODO
class MinimalityConnectedComponentChecker : public PartiallyDefinedGraphChecker
{
    vector<vertex_ordering_t> vertexOrderings;
    minimalit_check_config_t config;

    MinimalityConnectedComponentChecker(int frequency, partition_t initial_partition, vector<vertex_ordering_t> vertexOrderings, int cutoff)
    {
        this->name = "MinimalityChecker";
        this->frequency = frequency;
        this->config.initial_partition = initial_partition;
        this->vertexOrderings = vertexOrderings;
        this->config.cutoff = cutoff;
    }

    void checkProperty(const adjacency_matrix_t &matrix);

    if (checkConnectedComponents) // special check
    {
      // check swapping paritions with same size
      int N = intervallsColoring.size(); // last partition is only the single vertex
      for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
          if (intervallsColoring[i].second - intervallsColoring[i].first == intervallsColoring[j].second - intervallsColoring[j].first)
          {
            // check swapping this two blocks
            vertex_ordering_t vertex_ordering(vertices);
            for (int i = 0; i < vertices; i++)
            {
              vertex_ordering[i] = i;
            }

            for (int k = 0; k < intervallsColoring[i].second - intervallsColoring[i].first; k++)
            {
              swap(vertex_ordering[k + intervallsColoring[i].first], vertex_ordering[k + intervallsColoring[j].first]);
            }
            checkMinimality(matrix, vertex_ordering, {.initial_partition = initialPartition, .cutoff = cutoff});
          }

      // check swapping connected components
      checkSwaps(matrix);
    }
};
*/

#endif
