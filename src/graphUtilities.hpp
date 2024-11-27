#ifndef GRAPH_UTILITIES_CHECKER_H
#define GRAPH_UTILITIES_CHECHER_H

#include <vector>
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
#include <cstdint>
#include <map>
#include "useful.h"

using namespace boost;
using namespace std;

template <typename ConnectedMap>
class connected_recorder : public default_bfs_visitor
{
public:
    connected_recorder(ConnectedMap conn) : c(conn) {}

    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph& g) const
    {
        c[u] = true;
    }

private:
    ConnectedMap c;
};

template <typename ConnectedMap>
connected_recorder<ConnectedMap> record_connected(ConnectedMap d)
{
    return connected_recorder<ConnectedMap> (d);
}

class GraphUtilities 
{    
    public:
        GraphUtilities () : m_graph_size(2){}
        GraphUtilities (int graph_size) : m_graph_size(graph_size){}
    protected:
        int m_graph_size;
        boost_graph matrix2BoostGraph(const adjacency_matrix_t &matrix);
        void random_spanning_tree(const boost_graph &g, std::vector<Vertex> &predecessors);
        bool is_connected(const boost_graph &g);
        void shortest_path(const boost_graph &g, Vertex dest, std::vector<Vertex> &predecessors);
        vector<signed_edge_t> spanningTree2Clause(const vector<edge_t> &edges_spanning_tree);
};

#endif