#ifndef DECOMPOSABILITY_CHECKER_H
#define DECOMPOSABILITY_CHECKER_H
#include <vector>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <boost/graph/adjacency_list.hpp>
#include <boost/config.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/array.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <array>
#include <utility>
#include <random>
#include <ctime>
#include <cstdint>
#include <map>
#include "useful.h"
#include "graphChecker.hpp"
#include "graphUtilities.hpp"

using namespace boost;
using namespace std;

class SpanningTreeChecker : public PartiallyDefinedGraphChecker, public GraphUtilities
{
    public:
        SpanningTreeChecker(int graph_size)
        {
            name = "SpanningTreeChecker";
            m_graph_size = graph_size;
        }
        SpanningTreeChecker(int graph_size, int freq)
        {
            name = "SpanningTreeChecker";
            frequency = freq;
            m_graph_size = graph_size;
        }
        void checkProperty(const adjacency_matrix_t &matrix);
};

class ThreeDecomposabilityChecker : public FullyDefinedGraphChecker, public GraphUtilities
{
    private:
        bool m_full_search;
        int m_iteration_count_max = 10;

    public:
        ThreeDecomposabilityChecker (int graph_size, bool full_search)
        {
            name = "ThreeDecomposabilityChecker";
            m_graph_size = graph_size;
            m_full_search = full_search; 
        };

        ThreeDecomposabilityChecker (int graph_size, bool full_search, int iteration_count_max)
        {
            name = "ThreeDecomposabilityChecker";
            m_graph_size = graph_size;
            m_full_search = full_search; 
            m_iteration_count_max = iteration_count_max;
        };
        void checkProperty(const adjacency_matrix_t &matrix);

    protected:
        const static vector<vector<int>> possible_coloring_2;
        const static vector<vector<int>> possible_coloring_3;
        const static vector<vector<int>> possible_coloring_3_no_red;

        bool find_non_separating_arc(boost_graph &original_g, boost_graph &g, Vertex s, Vertex t, vector<int> &non_separating_cycle);
        bool find_non_separating_cycle(boost_graph &g);
        vector<edge_t> preprocess_graph(boost_graph &g);
        void move_edge_forward(Vertex v_2, Vertex v_1, boost_graph &g, boost_graph &spanning_tree, boost_graph &tree_complement);
        std::pair<bool, std::pair<Vertex, Vertex>> find_next_edge(const boost_graph &g);
        vector<edge_t> three_decomposition_heuristic(boost_graph g);
        bool coloring_compatible(const vector<int>& old_coloring, const vector<int>& new_coloring, int length);
        void combine_components(Vertex smaller, Vertex bigger, map<Vertex, Vertex> &vertex_comp_dict, map<Vertex, Vertex> &vertex_comp_change);
        bool try_compatible_coloring(VertexIt vit, const boost_graph &g, map<Edge, int> &edge_color_dict, map<Vertex, Vertex> &vertex_comp_dict, const vector<vector<int>> &possible_coloring);
        bool assign_edge_coloring(VertexIt vit, const boost_graph &g, map<Edge, int> &edge_color_dict, map<Vertex, Vertex> &vertex_comp_dict);
        vector<edge_t> three_decomposition_full(boost_graph &g);
};          

#endif