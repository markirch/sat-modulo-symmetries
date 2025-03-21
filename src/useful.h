#ifndef USEFUL_H
#define USEFUL_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <string>

using std::vector;
// using std::pair;

typedef int lit_t;
typedef vector<lit_t> clause_t;
typedef vector<clause_t> cnf_t;
typedef vector<vector<lit_t>> edge_vars_t;

enum LogLevel
{
    LOG_LEVEL_INFO,
    LOG_LEVEL_WARNING,
    LOG_LEVEL_ERROR,
    LOG_LEVEL_DEBUG
};

#define LOG_LEVEL LOG_LEVEL_INFO

#define LOG(level, message)                  \
    if (level <= LOG_LEVEL)            \
    {                                        \
        std::ostringstream oss;              \
        oss << message;                      \
        std::cout << oss.str() << std::endl; \
    }

#define PRINT_CURRENT_LINE                                \
    if (LOG_LEVEL >= LOG_LEVEL_DEBUG)               \
    {                                                     \
        printf("Line %d, file %s, function %s\n", __LINE__, __FILE__, __func__); \
        fflush(stdout);                                   \
    }

#define EXIT_UNWANTED_STATE                                                          \
    {                                                                                \
        printf("Error: unexpected state at line %d, file %s\n", __LINE__, __FILE__); \
        exit(EXIT_FAILURE);                                                          \
    }

typedef int vertex_t;
typedef std::pair<vertex_t, vertex_t> edge_t;

typedef enum
{
    truth_value_unknown = 2,
    truth_value_true = 1,
    truth_value_false = 0
} truth_value_t;

typedef std::vector<std::vector<truth_value_t>> adjacency_matrix_t;
typedef std::vector<std::vector<truth_value_t>> incidence_matrix_t;
typedef std::pair<truth_value_t, edge_t> signed_edge_t;

typedef vector<signed_edge_t> forbidden_graph_t; // forbidden graph as signed edge list.

void printAdjacencyMatrix(const adjacency_matrix_t &matrix, bool printFullMatrix);
void printPartiallyDefinedAdjacencyMatrix(const adjacency_matrix_t &matrix);
void printHypergraph(const adjacency_matrix_t &matrix, int nVertices[2]);

adjacency_matrix_t getIntersectionMatrix(const adjacency_matrix_t &incidence_matrix, int b_vertices[2]);

void file2cnf(std::ifstream &file, cnf_t &cnf, int &maxVar);

#endif
