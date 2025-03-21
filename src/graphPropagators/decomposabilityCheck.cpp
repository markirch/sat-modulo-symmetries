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
#include "decomposabilityCheck.hpp"

void print_edges(const vector<edge_t> &edge_list)
{
    auto it = edge_list.begin();
    cout << '[';
    while(it != edge_list.end())
    {
        cout << '(' << it->first << ',' << it->second << ')';
        if (it != edge_list.end() - 1)
            cout << ", ";
        ++it;
    }
    cout << ']' << endl;
}

void SpanningTreeChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    boost_graph g = matrix2BoostGraph(matrix);
    if (!is_connected(g))
       return;

    //std::vector<int> predecessors(m_graph_size, -1);
    std::vector<Vertex> predecessors(m_graph_size);
    random_spanning_tree(g, predecessors);

    vector<edge_t> edges_spanning_tree;
    for (int i = 1; i < m_graph_size; ++i)
    {
        // if (predecessors[i] == -1)
        //     return;
        edges_spanning_tree.push_back(std::make_pair(i, predecessors[i]));
    }

    //cout << calls << ',' << frequency << endl;
    //print_edges(edges_spanning_tree);
    auto forbidden_tree = spanningTree2Clause(edges_spanning_tree);
    for (auto it = forbidden_tree.begin(); it != forbidden_tree.end(); ++it)
    {
        if (it->first == truth_value_true)
            continue;
        auto i = (it->second).first, j = (it->second).second;
        if (matrix[i][j] == truth_value_true)
            return;
    }
    throw forbidden_tree;
}

const vector<vector<int>> ThreeDecomposabilityChecker::possible_coloring_2 = {{edge_color_green, edge_color_blue},
                                                                            {edge_color_blue, edge_color_green},
                                                                            {edge_color_green, edge_color_green}};
const vector<vector<int>> ThreeDecomposabilityChecker::possible_coloring_3 = {{edge_color_red, edge_color_red, edge_color_green}, 
                                                                            {edge_color_red, edge_color_green, edge_color_red},
                                                                            {edge_color_green, edge_color_red, edge_color_red},  
                                                                            {edge_color_green, edge_color_blue, edge_color_green},
                                                                            {edge_color_green, edge_color_green, edge_color_blue},
                                                                            {edge_color_blue, edge_color_green, edge_color_green},
                                                                            {edge_color_green, edge_color_green, edge_color_green}};
const vector<vector<int>> ThreeDecomposabilityChecker::possible_coloring_3_no_red = {{edge_color_green, edge_color_blue, edge_color_green},
                                                                                   {edge_color_green, edge_color_green, edge_color_blue},
                                                                                   {edge_color_blue, edge_color_green, edge_color_green},
                                                                                   {edge_color_green, edge_color_green, edge_color_green}}; 

bool ThreeDecomposabilityChecker::find_non_separating_arc(boost_graph &original_g, boost_graph &g, Vertex s, Vertex t, vector<int> &non_separating_cycle)
{
    if (s == t)
    {
        if(is_connected(original_g))
            return true;
        return false;
    }
    graph_traits<boost_graph>::out_edge_iterator eit, eend;
    std::tie(eit, eend) = boost::out_edges(s, g);
    vector<Vertex> next_candidates;
    while(eit != eend)
    {
        Vertex s_new = boost::target(*eit, g);
        if (non_separating_cycle[s_new] == -1)
            next_candidates.push_back(s_new);
        ++eit;
    }
    auto it = next_candidates.begin();
    while(it != next_candidates.end())
    {
        Vertex s_new = *it;
        non_separating_cycle[s_new] = s;
        boost::remove_edge(s, s_new, g);
        boost::remove_edge(s, s_new, original_g);
        if (find_non_separating_arc(original_g, g, s_new, t, non_separating_cycle))
            return true;
        boost::add_edge(s, s_new, g);
        boost::add_edge(s, s_new, original_g);
        non_separating_cycle[s_new] = -1;
        ++it;
    }
    return false;
}

bool ThreeDecomposabilityChecker::find_non_separating_cycle(boost_graph &g)
{
    
    // Do heuristic search here
    vector<Vertex> edges_spanning_tree(m_graph_size, Vertex());
    random_spanning_tree(g, edges_spanning_tree);
    boost_graph spanning_tree, g_copy(g);
    
    VertexIt vit, vend;
    std::tie(vit, vend) = boost::vertices(g);
    vit++;
    while(vit != vend)
    {
        boost::add_edge(*vit, edges_spanning_tree[*vit], spanning_tree);
        boost::remove_edge(*vit, edges_spanning_tree[*vit], g_copy);
        ++vit;
    }
    
    graph_traits<boost_graph>::edge_iterator eit, eend;
    std::tie(eit, eend) = boost::edges(g_copy);
    while(eit != eend)
    { 
        Vertex s = boost::source(*eit, g_copy), t = boost::target(*eit, g_copy);
        vector<Vertex> predecessor(m_graph_size);
        shortest_path(spanning_tree, s, predecessor);
        Vertex t_new = t;
        boost::remove_edge(s, t, g);
        while(true)
        {
            boost::remove_edge(t_new, predecessor[t_new], g);
            t_new = predecessor[t_new];
            if (t_new == s)
                break;
        }
        if (is_connected(g))
        
            return true;
        t_new = t;
        while(true)
        {
            boost::add_edge(t_new, predecessor[t_new], g);
            t_new = predecessor[t_new];
            if (t_new == s)
                break;
        }
        boost::add_edge(s, t, g);
        ++eit;
    }

    // Do full search here
    if (m_full_search)
    {
        graph_traits<boost_graph>::edge_iterator eit, eend;
        boost_graph g_copy(g);
        std::tie(eit, eend) = boost::edges(g_copy);
        while(eit != eend)
        {
            vector<int> non_separating_cycle(boost::num_vertices(g), -1);
            Vertex s = boost::source(*eit, g_copy), t = boost::target(*eit, g_copy);
            non_separating_cycle[s] = t;
            boost::remove_edge(s, t, g_copy);
            boost::remove_edge(s, t, g);

            if (find_non_separating_arc(g, g_copy, s, t, non_separating_cycle))
                return true;
            boost::add_edge(s, t, g);
            std::tie(eit, eend) = boost::edges(g_copy);
        }
    }
    return false;
}

vector<edge_t> ThreeDecomposabilityChecker::preprocess_graph(boost_graph &g)
{
    vector<edge_t> safe_edges;
    VertexIt vit, vend;
    graph_traits<boost_graph>::adjacency_iterator ait, aend;
    vector<Vertex> deg_1_list;
    vertex_name_map_t vertex_name_map_g = get(vertex_name, g);

    while(find_non_separating_cycle(g))
    {
        std::tie(vit, vend) = boost::vertices(g);
        while(vit != vend)
        {
            if (boost::out_degree(*vit, g) == 1)
            {
                std::tie(ait, aend) = boost::adjacent_vertices(*vit, g);
                safe_edges.push_back(make_pair(vertex_name_map_g[*vit], vertex_name_map_g[*ait]));
                boost::clear_vertex(*vit, g);
                boost::remove_vertex(*vit, g);
                std::tie(vit, vend) = boost::vertices(g);
                continue;
            }
            ++vit;
        }
    }
    return safe_edges;
}

void ThreeDecomposabilityChecker::move_edge_forward(Vertex v_2, Vertex v_1, boost_graph &g, boost_graph &spanning_tree, boost_graph &tree_complement)
{
    boost::remove_edge(v_2, v_1, tree_complement);
    vector<Vertex> predecessor(m_graph_size, -1);
    shortest_path(spanning_tree, v_2, predecessor);
    Vertex v_new = v_1;
    while(true)
    {
        boost::remove_edge(v_new, predecessor[v_new], g);
        v_new = predecessor[v_new];
        if (v_new == v_2)
            break;
    }
    boost::remove_edge(v_2, v_1, g);
    if (is_connected(g))
    {
        vector<Vertex> edges_spanning_tree(m_graph_size, -1);
        random_spanning_tree(g, edges_spanning_tree);
        tree_complement = g;
        spanning_tree.clear();
        VertexIt vit, vend;
        std::tie(vit, vend) = boost::vertices(g);
        vit++;
        while(vit != vend)
        {
            boost::add_edge(*vit, edges_spanning_tree[*vit], spanning_tree);
            boost::remove_edge(*vit, edges_spanning_tree[*vit], tree_complement);
            ++vit;
        }
        return;
    }

    boost::add_edge(v_2, v_1, g);
    boost::add_edge(v_2, v_1, spanning_tree);
    v_new = v_1;
    while(true)
    {
        boost::add_edge(v_new, predecessor[v_new], g);
        v_new = predecessor[v_new];
        if (v_new == v_2)
            break;
    }
    v_new = v_1;
    while(true)
    {
        if (boost::out_degree(v_new, tree_complement) == 0 && boost::out_degree(predecessor[v_new], tree_complement) == 0)
        {
            boost::remove_edge(v_new, predecessor[v_new], spanning_tree);
            boost::add_edge(v_new, predecessor[v_new], tree_complement);
            return;
        }
        v_new = predecessor[v_new];
        if (v_new == v_2)
            break;
    }

    boost::remove_edge(v_1, predecessor[v_1], spanning_tree);
    boost::add_edge(v_1, predecessor[v_1], tree_complement);
    return;
}

std::pair<bool, std::pair<Vertex, Vertex>> ThreeDecomposabilityChecker::find_next_edge(const boost_graph &g)
{
    graph_traits<boost_graph>::edge_iterator eit, eend;
    std::tie(eit, eend) = boost::edges(g);
    Vertex s, t;
    while (eit != eend)
    {
        s = boost::source(*eit, g);
        t = boost::target(*eit, g);
        if (boost::out_degree(s, g) == 2 && boost::out_degree(t, g) == 1)
            return make_pair(true, make_pair(s, t));
        if (boost::out_degree(s, g) == 1 && boost::out_degree(t, g) == 2)
            return make_pair(true, make_pair(t, s));
        ++eit;
    }
    return make_pair(false, make_pair(Vertex(), Vertex()));
}

vector<edge_t> ThreeDecomposabilityChecker::three_decomposition_heuristic(boost_graph g)
{
    vector<Vertex> edges_spanning_tree(m_graph_size);
    random_spanning_tree(g, edges_spanning_tree);
    boost_graph tree_complement(g);
    boost_graph spanning_tree;
    VertexIt vit, vend;
    std::tie(vit, vend) = boost::vertices(g);
    vit++;
    while(vit != vend)
    {
        boost::add_edge(*vit, edges_spanning_tree[*vit], spanning_tree);
        boost::remove_edge(*vit, edges_spanning_tree[*vit], tree_complement);
        ++vit;
    }

    bool found;
    Vertex v_2, v_1;
    std::pair<Vertex, Vertex> v_pair;

    int iteration_count = 0;
    std::tie(found, v_pair) = find_next_edge(tree_complement);
    std::tie(v_2, v_1) = v_pair;
    while (found && iteration_count < m_iteration_count_max)
    {
        iteration_count++;
        move_edge_forward(v_2, v_1, g, spanning_tree, tree_complement);
        std::tie(found, v_pair) = find_next_edge(tree_complement);
        std::tie(v_2, v_1) = v_pair;
    }
    if (found)
        return vector<edge_t>();
    vector<edge_t> final_edges_spanning_tree;
    graph_traits<boost_graph>::edge_iterator eit, eend;
    std::tie(eit, eend) = boost::edges(spanning_tree);
    vertex_name_map_t vertex_name_map_g = get(vertex_name, g);
    while(eit != eend)
    {
        final_edges_spanning_tree.push_back(make_pair(vertex_name_map_g[boost::source(*eit, spanning_tree)], vertex_name_map_g[boost::target(*eit, spanning_tree)]));
        ++eit;
    }
    return final_edges_spanning_tree;
}

bool ThreeDecomposabilityChecker::coloring_compatible(const vector<int>& old_coloring, const vector<int>& new_coloring, int length)
{
    for(int i = 0; i < length; ++i)
    {
        if (old_coloring[i] != 0 && old_coloring[i] != new_coloring[i])
            return false;
    }
    return true;
}

void ThreeDecomposabilityChecker::combine_components(Vertex smaller, Vertex bigger, map<Vertex, Vertex> &vertex_comp_dict, map<Vertex, Vertex> &vertex_comp_change)
{
    if (smaller > bigger)
    {
        combine_components(bigger, smaller, vertex_comp_dict, vertex_comp_change);
        return;
    }

    for (auto it = vertex_comp_dict.begin(); it != vertex_comp_dict.end(); ++it)    
    {
        if (it->second == bigger)
        {
            if (vertex_comp_change.find(it->first) == vertex_comp_change.end())
                vertex_comp_change[it->first] = it->second;
            vertex_comp_dict[it->first] = smaller;
        }
    }
}

bool ThreeDecomposabilityChecker::try_compatible_coloring(VertexIt vit, const boost_graph &g, map<Edge, int> &edge_color_dict, map<Vertex, Vertex> &vertex_comp_dict, const vector<vector<int>> &possible_coloring)
{
    vector<Vertex> neighbors;
    vector<int> neighbor_coloring;
    int num_neighbors = 0;
    graph_traits<boost_graph>::adjacency_iterator ait, aend;
    std::tie(ait, aend) = boost::adjacent_vertices(*vit, g);
    //cout << "Neighbor Info: ";
    while(ait != aend)
    {
        neighbors.push_back(*ait);
        neighbor_coloring.push_back(edge_color_dict[make_pair(*vit, *ait)]);
        ++ait;
        ++num_neighbors;
    }
    //cout << endl;
    for (auto it = possible_coloring.begin(); it != possible_coloring.end(); ++it)
    {
        map<Vertex, Vertex> vertex_comp_change;
        vector<Edge> edge_color_change;
        if(!coloring_compatible(neighbor_coloring, *it, num_neighbors))
            continue;
        bool green_loop = false;
        for(int i = 0; i < num_neighbors; ++i)
        {
            if (neighbor_coloring[i] == unassigned)
            {
                if ((*it)[i] == edge_color_green)
                {
                    if(vertex_comp_dict[*vit] == vertex_comp_dict[neighbors[i]])
                    {
                        green_loop = true;
                        break;
                    }
                    combine_components(vertex_comp_dict[*vit], vertex_comp_dict[neighbors[i]], vertex_comp_dict, vertex_comp_change);
                }
                edge_color_change.push_back(make_pair(*vit, neighbors[i]));
                edge_color_dict[make_pair(*vit, neighbors[i])] = (*it)[i];
                edge_color_dict[make_pair(neighbors[i],*vit)] = (*it)[i];
            }
        }
        if (green_loop)
        {
            for(auto vit = vertex_comp_change.begin(); vit != vertex_comp_change.end(); ++vit)
                vertex_comp_dict[vit->first] = vit->second;
            for(auto eit = edge_color_change.begin(); eit != edge_color_change.end(); ++eit)
            {
                edge_color_dict[*eit] = 0;
                edge_color_dict[make_pair(eit->second, eit->first)] = 0;
            }
            continue;
        }

        if(assign_edge_coloring(vit + 1, g, edge_color_dict, vertex_comp_dict))
            return true;
        for(auto vit = vertex_comp_change.begin(); vit != vertex_comp_change.end(); ++vit)
            vertex_comp_dict[vit->first] = vit->second;
        for(auto eit = edge_color_change.begin(); eit != edge_color_change.end(); ++eit)
        {
            edge_color_dict[*eit] = 0;
            edge_color_dict[make_pair(eit->second, eit->first)] = 0;
        }    
    }
    return false;
}

bool ThreeDecomposabilityChecker::assign_edge_coloring(VertexIt vit, const boost_graph &g, map<Edge, int> &edge_color_dict, map<Vertex, Vertex> &vertex_comp_dict)
{
    VertexIt it, vend;
    std::tie(it, vend) = boost::vertices(g);
    VertexIt vbeg = it;
    if (vit == vend)
    {
        while(it != vend)
        {
            if(vertex_comp_dict[*it] != *vbeg)
                return false;
            ++it;
        }
        return true; 
    }

    if (boost::out_degree(*vit, g) == 2)
        return try_compatible_coloring(vit, g, edge_color_dict, vertex_comp_dict, possible_coloring_2);
    if (m_full_search)
        return try_compatible_coloring(vit, g, edge_color_dict, vertex_comp_dict, possible_coloring_3_no_red);
    return try_compatible_coloring(vit, g, edge_color_dict, vertex_comp_dict, possible_coloring_3);
}

vector<edge_t> ThreeDecomposabilityChecker::three_decomposition_full(boost_graph &g)
{
    map<Edge, int> edge_color_dict;
    map<Vertex, Vertex> vertex_comp_dict;
    VertexIt vit, vend;
    std::tie(vit, vend) = boost::vertices(g);
    while(vit != vend)
    {
        vertex_comp_dict[*vit] = *vit;
        ++vit;
    }
    graph_traits<boost_graph>::edge_iterator eit, eend;
    std::tie(eit, eend) = boost::edges(g);
    while (eit != eend)
    {
        edge_color_dict[make_pair(boost::source(*eit, g), boost::target(*eit, g))] = unassigned;
        edge_color_dict[make_pair(boost::target(*eit, g), boost::source(*eit, g))] = unassigned;
        ++eit;
    }
    std::tie(vit, vend) = boost::vertices(g);
    if (assign_edge_coloring(vit, g, edge_color_dict, vertex_comp_dict))
    {
        vector<edge_t> edges_spanning_tree;
        vertex_name_map_t vertex_name_map_g = get(vertex_name, g);
        for (auto it = edge_color_dict.begin(); it != edge_color_dict.end(); ++it)
        {
            if((it->first).first > (it->first).second)
                continue;
            if (it->second != edge_color_green)
                continue;
            edges_spanning_tree.push_back(make_pair(vertex_name_map_g[(it->first).first], vertex_name_map_g[(it->first).second]));
        }
        return edges_spanning_tree;
    }
    return vector<edge_t>();  
}

void ThreeDecomposabilityChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    boost_graph g = matrix2BoostGraph(matrix);
    vector<edge_t> safe_edges = preprocess_graph(g);
    if (boost::num_vertices(g) > 1)
    {
        vector<edge_t> spanning_tree = three_decomposition_heuristic(g);
        if (spanning_tree.size() < 1)
            spanning_tree = three_decomposition_full(g);
        if (spanning_tree.size() < 1)
            return;
        safe_edges.insert(safe_edges.end(), spanning_tree.begin(), spanning_tree.end());
    }
    //print_edges(safe_edges);
    throw spanningTree2Clause(safe_edges);
}