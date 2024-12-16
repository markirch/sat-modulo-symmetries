#include <cstdlib>
#include <iostream>
#include <ctime>
#include <utility>
#include <algorithm>
#include <fstream>
#include <queue>
#include <array>
#include <utility>
#include <random>
#include <ctime>
#include <cstdint>
#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/config.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/array.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
//#include "useful.h"
#include "graphUtilities.hpp"

boost_graph GraphUtilities::matrix2BoostGraph(const adjacency_matrix_t &matrix)
{
    boost_graph g;
    vertex_name_map_t vertex_name_map_g = get(vertex_name, g);
    int i, j;
    for (i = 0; i < m_graph_size; ++i)
    {
        Vertex u = add_vertex(g);
        vertex_name_map_g[u] = u;
    }
    VertexIt vit, vend, vitp;
    std::tie(vit, vend) = boost::vertices(g);
    i = 0;
    while(vit != vend)
    {
        j = i + 1;
        vitp = vit + 1;
        while(vitp != vend)
        {
            if (matrix[i][j] == truth_value_true)
                boost::add_edge(*vit, *vitp, g);
            ++j;
            ++vitp;
        }
        ++i;
        ++vit;
    }
    return g;
}

void GraphUtilities::random_spanning_tree(const boost_graph &g, std::vector<Vertex> &predecessors)
{
    std::mt19937 gen{static_cast<uint32_t>(std::time(0))};
    Vertex src = Vertex();
    boost::random_spanning_tree(g, gen, boost::predecessor_map(boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, g))).root_vertex(src));
}

bool GraphUtilities::is_connected(const boost_graph &g)
{
    Vertex src = 0;
    auto nof_nodes = boost::num_vertices(g);
    std::vector<bool> connected(nof_nodes, false);
    boost::breadth_first_search(g, src, boost::visitor(record_connected(boost::make_iterator_property_map(connected.begin(), boost::get(boost::vertex_index, g)))));
    for (int i = 1; i < nof_nodes; ++i)
        if (!connected[i])
            return false;
    return true;
}

void GraphUtilities::shortest_path(const boost_graph &g, Vertex dest, std::vector<Vertex> &predecessors)
{
    boost::breadth_first_search(g, dest, 
        boost::visitor(boost::make_bfs_visitor(boost::record_predecessors(boost::make_iterator_property_map(predecessors.begin(), boost::get(boost::vertex_index, g)), boost::on_tree_edge{}))));
}

vector<signed_edge_t> GraphUtilities::spanningTree2Clause(const vector<edge_t> &edges_spanning_tree)
{
    boost_graph spanning_tree(edges_spanning_tree.begin(), edges_spanning_tree.end(), m_graph_size);
    vertex_name_map_t vertex_name_map_g = get(vertex_name, spanning_tree);
    VertexIt vit, vend;
    std::tie(vit, vend) = boost::vertices(spanning_tree);
    while (vit != vend)
    {
        vertex_name_map_g[*vit] = *vit;
        ++vit;
    }

    std::tie(vit, vend) = boost::vertices(spanning_tree);
    vector<Vertex> degree_1_list, degree_2_list;
    while(vit != vend)
    {
        boost::graph_traits<boost_graph>::degree_size_type degree_type = out_degree(*vit, spanning_tree);
        if (degree_type == 1)
        {
            degree_1_list.push_back(*vit);
            ++vit;
            continue;
        }
        if (degree_type == 2)
            degree_2_list.push_back(*vit);
        ++vit;
    }

    vector<edge_t> error_pairs;
    vector<Vertex>::iterator it_1 = degree_1_list.begin();
    while (it_1 != degree_1_list.end())
    {
        vector<Vertex>::iterator it_2 = degree_2_list.begin();
        while (it_2 != degree_2_list.end())
        {
            if (find(edges_spanning_tree.begin(), edges_spanning_tree.end(), make_pair(vertex_name_map_g[*it_1], vertex_name_map_g[*it_2])) == edges_spanning_tree.end() 
            && find(edges_spanning_tree.begin(), edges_spanning_tree.end(), make_pair(vertex_name_map_g[*it_2], vertex_name_map_g[*it_1])) == edges_spanning_tree.end())
                error_pairs.push_back(make_pair(vertex_name_map_g[*it_1], vertex_name_map_g[*it_2]));
            ++ it_2;
        }
        ++ it_1;
    }
    vector<signed_edge_t> clause;
    for (auto it = edges_spanning_tree.begin(); it != edges_spanning_tree.end(); ++it)
        clause.push_back(make_pair(truth_value_true, *it));
    for (auto it = error_pairs.begin(); it != error_pairs.end(); ++it)
        clause.push_back(make_pair(truth_value_false, *it));
    return clause;
}
