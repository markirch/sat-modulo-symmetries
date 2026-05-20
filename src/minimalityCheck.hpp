#ifndef MIN_CHECK_H
#define MIN_CHECK_H

#include "useful.h"
#include "graphChecker.hpp"

#define DEFAULT_MINIMALITY_FREQUENCY 30
#define DEFAULT_MINIMALITY_CUTOFF 200000

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

struct minimality_config_t
{
    int frequency = DEFAULT_MINIMALITY_FREQUENCY;
    partition_t initialPartition;
    int cutoff = DEFAULT_MINIMALITY_CUTOFF;
    vector<vertex_ordering_t> vertexOrderings;
    bool turnoffSMS = false;
    bool useColexOrdering = false;
    std::string symBreakClausesFile;
};

inline partition_t getDefaultInitialPartition(int vertices)
{
    partition_t partition(vertices, false);
    partition[0] = true;
    return partition;
}

inline vector<vertex_ordering_t> getDefaultVertexOrderings(int vertices)
{
    vector<vertex_ordering_t> vertexOrderings;
    vertex_ordering_t vertexOrdering(vertices, 0);
    for (int i = 0; i < vertices; i++)
        vertexOrdering[i] = i;
    vertexOrderings.push_back(vertexOrdering);
    return vertexOrderings;
}

/**
 * Check if graph is minimal and add a clause if this is not the case.
 * Throws std::vector<signed_edge_t> if current partially defined graph is not minimal
 * Throws LimitReachedException if cutoff limit is reached
 */
void checkMinimality(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityColex(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityDir(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityComplement(adjacency_matrix_t &adjacency_matrix, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);
void checkMinimalityMultiple(vector<adjacency_matrix_t> &adjacency_matrices, vertex_ordering_t vertex_ordering, minimalit_check_config_t config);

// basic minimality check with potentially different initial vertexOrderings
class MinimalityChecker : public PartiallyDefinedGraphChecker
{
private:
    vector<vertex_ordering_t> vertexOrderings;
    minimalit_check_config_t config;
    FILE *symBreakClauses = nullptr;
    bool ownsSymBreakClauses = false;

    bool directed = false;
    bool maximize = false;
    bool useColexOrdering = false;

public:
    MinimalityChecker(int frequency, partition_t initial_partition, vector<vertex_ordering_t> vertexOrderings, int cutoff, FILE *symBreakClauses)
    {
        this->name = "MinimalityChecker";
        this->frequency = frequency;
        this->config.initial_partition = initial_partition;
        this->vertexOrderings = vertexOrderings;
        this->config.cutoff = cutoff;
        this->symBreakClauses = symBreakClauses;
        this->useColexOrdering = false;
    }

    MinimalityChecker(struct minimality_config_t config, int vertices, bool directed = false, bool maximize = false)
    {
        this->name = "MinimalityChecker";
        this->frequency = config.frequency;
        this->config.initial_partition = config.initialPartition;
        this->vertexOrderings = config.vertexOrderings;
        this->config.cutoff = config.cutoff;
        this->useColexOrdering = config.useColexOrdering;
        if (!config.symBreakClausesFile.empty())
        {
            this->symBreakClauses = fopen(config.symBreakClausesFile.c_str(), "w");
            if (!this->symBreakClauses)
                EXIT_UNWANTED_STATE
            this->ownsSymBreakClauses = true;
            fprintf(this->symBreakClauses, "{ \"sym_clauses\": [\n");
        }

        // handle uninitialized values
        if (this->config.initial_partition.empty())
            this->config.initial_partition = getDefaultInitialPartition(vertices);

        if ((int)this->config.initial_partition.size() != vertices)
            EXIT_UNWANTED_STATE

        if (this->vertexOrderings.empty())
            this->vertexOrderings = getDefaultVertexOrderings(vertices);

        this->directed = directed;
        this->maximize = maximize;
    }

    ~MinimalityChecker()
    {
        if (ownsSymBreakClauses && symBreakClauses)
        {
            fprintf(symBreakClauses, " ]\n}\n");
            fclose(symBreakClauses);
        }
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

#endif
