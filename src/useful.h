#ifndef USEFUL_H
#define USEFUL_H

#include <vector>
#include <fstream>
#include <algorithm>
#include <cstdint>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>

/*
#include <algorithm>
#include <utility>
#include <iostream>
#include <numeric>
#include <string>
#include <string.h>
#include <cassert>
*/

using std::vector;
//using std::pair;

typedef int lit_t;
typedef vector<lit_t> clause_t;
typedef vector<clause_t> cnf_t;
typedef vector<vector<lit_t>> edge_vars_t;

#define PRINT_CURRENT_LINE                            \
    printf("Line %d, file %s\n", __LINE__, __FILE__); \
    fflush(stdout);

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

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_name_t, vertex_t>> boost_graph;
typedef boost::graph_traits<boost_graph>::vertex_descriptor Vertex;
typedef std::pair<Vertex, Vertex> Edge;
typedef boost::property_map<boost_graph, boost::vertex_name_t>::type vertex_name_map_t;
typedef boost::graph_traits<boost_graph>::vertex_iterator VertexIt;

enum{unassigned, edge_color_green, edge_color_red, edge_color_blue};

void printAdjacencyMatrix(const adjacency_matrix_t &matrix, bool printFullMatrix);
void printPartiallyDefinedAdjacencyMatrix(const adjacency_matrix_t &matrix);
void printHypergraph(const adjacency_matrix_t &matrix, int nVertices[2]);

adjacency_matrix_t getIntersectionMatrix(const adjacency_matrix_t &incidence_matrix, int b_vertices[2]);

void file2cnf(std::ifstream &file, cnf_t &cnf, int &maxVar);

#endif
