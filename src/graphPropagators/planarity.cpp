#include <boost/graph/boyer_myrvold_planar_test.hpp>
//#include <boost/graph/is_kuratowski_subgraph.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "planarity.hpp"

using std::pair;
using std::make_pair;
using std::min, std::max;
//using std::cout, std::endl;

vector<pair<int, int>> deleteVerticesWithDegree1(vector<pair<int, int>> &edges)
{
    std::unordered_map<int, int> degreeCount;

    // Count the degree of each vertex
    for (const auto &edge : edges)
    {
        degreeCount[edge.first]++;
        degreeCount[edge.second]++;
    }

    // Create a new edge list without vertices of degree 1
    vector<pair<int, int>> newEdges;
    for (const auto &edge : edges)
    {
        if (degreeCount[edge.first] > 1 && degreeCount[edge.second] > 1)
        {
            newEdges.push_back(edge);
        }
        else
        {
            degreeCount[edge.first]--;
            degreeCount[edge.second]--;
        }
    }

    return newEdges;
}

vector<pair<int, int>> deleteVerticesWithDegree1All(vector<pair<int, int>> edges)
{
    auto edgesNew = deleteVerticesWithDegree1(edges);
    while (edgesNew.size() < edges.size())
    {
        edges = edgesNew;
        edgesNew = deleteVerticesWithDegree1(edges);
    }
    return edgesNew;
}

// returns a K_5 or K_{3,3} subgraph if not planar, otherwise empty list
vector<pair<int, int>> testPlanarity(const adjacency_matrix_t &m)
{

    int n = (int)m.size(); // vertices is already used in the boost namespace
    using namespace boost;
    typedef adjacency_list<vecS,
                           vecS,
                           undirectedS,
                           property<vertex_index_t, int>,
                           property<edge_index_t, int>>
        graph;

    // Create graph in format suitable for boost
    graph g(n);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (m[i][j] == truth_value_true)
                add_edge(i, j, g);

    // Initialize the interior edge index
    property_map<graph, edge_index_t>::type e_index = get(edge_index, g);
    graph_traits<graph>::edges_size_type edge_count = 0;
    graph_traits<graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
        put(e_index, *ei, edge_count++);

    // Test for planarity
    typedef std::vector<graph_traits<graph>::edge_descriptor>
        kuratowski_edges_t;
    kuratowski_edges_t kuratowski_edges;
    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                                     boyer_myrvold_params::kuratowski_subgraph =
                                         std::back_inserter(kuratowski_edges)))
    {
        // std::cout << "Input graph is planar" << std::endl;
        return vector<pair<int, int>>();
    }
    else
    {
        // std::cout << "Input graph is not planar" << std::endl;
        // std::cout << "Edges in the Kuratowski subgraph: ";
        vector<pair<int, int>> kurEdges;
        kuratowski_edges_t::iterator ki, ki_end;
        ki_end = kuratowski_edges.end();
        for (ki = kuratowski_edges.begin(); ki != ki_end; ++ki)
        {
            // std::cout << *ki << " ";
            // printf(" %d %d", source(*ki, g), target(*ki, g));
            kurEdges.push_back(make_pair(source(*ki, g), target(*ki, g)));
        }
        // std::cout << std::endl;
        return deleteVerticesWithDegree1All(kurEdges);
    }
}

// ensure that graph is planar
void PlanarityChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    // printf("Check planarity\n");
    auto edges = testPlanarity(matrix);
    if (!edges.empty())
    {
        forbidden_graph_t forbiddenGraph;
        for (auto e : edges)
            forbiddenGraph.push_back(make_pair(truth_value_true, e));
        throw forbiddenGraph;
    }
}

// ensure that underlying graph is planar
void DirectedPlanarityChecker::checkProperty(const adjacency_matrix_t &m)
{
    auto underlyingMatrix = m;
    int vertices = (int)m.size();

    for (int i = 0; i < vertices; i++)
        for (int j = i + 1; j < vertices; j++)
        {
            if (m[i][j] == truth_value_true || m[j][i] == truth_value_true)
                underlyingMatrix[j][i] = underlyingMatrix[i][j] = truth_value_true;
        }

    auto edges = testPlanarity(underlyingMatrix);
    if (!edges.empty())
    {
        forbidden_graph_t forbiddenGraph;
        for (auto e : edges)
        {
            if (m[e.first][e.second] == truth_value_true)
                forbiddenGraph.push_back(make_pair(truth_value_true, e));
            else // edge in other direction present
                forbiddenGraph.push_back(make_pair(truth_value_true, make_pair(e.second, e.first)));
        }
        throw forbiddenGraph;
    }
}

// ensure that thickness is two where the decomposition is represented by a directed graph where edges in both direction belong to the first graph and all single edges to the other one. By minimalty one can assume that unidirectional are from higher to lower vertices.
void ThicknessTwoChecker::checkProperty(const adjacency_matrix_t &m)
{
    int vertices = (int)m.size();
    auto m1 = vector<vector<truth_value_t>>(vertices, vector<truth_value_t>(vertices, truth_value_false));
    for (int i = 0; i < vertices; i++)
        for (int j = i + 1; j < vertices; j++)
            if (m[i][j] == truth_value_true && m[j][i] == truth_value_true)
                m1[i][j] = m1[j][i] = truth_value_true;

    auto edges = testPlanarity(m1);
    if (!edges.empty())
    {
        forbidden_graph_t forbiddenGraph;
        for (auto e : edges)
        {
            forbiddenGraph.push_back(make_pair(truth_value_true, e));
            forbiddenGraph.push_back(make_pair(truth_value_true, make_pair(e.second, e.first))); // in both directions
        }
        throw forbiddenGraph;
    }

    // test edges only in lower triangular matrix
    auto m2 = vector<vector<truth_value_t>>(vertices, vector<truth_value_t>(vertices, truth_value_false));
    for (int i = 0; i < vertices; i++)
        for (int j = i + 1; j < vertices; j++)
            if (m[j][i] == truth_value_true && m[i][j] == truth_value_false) // only adjacent on lower side
                m2[i][j] = m2[j][i] = truth_value_true;
    edges = testPlanarity(m2);
    if (!edges.empty())
    {
        forbidden_graph_t forbiddenGraph;
        for (auto e : edges)
        {
            forbiddenGraph.push_back(make_pair(truth_value_true, make_pair(max(e.first, e.second), min(e.first, e.second))));  // from higher to lower vertex present
            forbiddenGraph.push_back(make_pair(truth_value_false, make_pair(min(e.first, e.second), max(e.first, e.second)))); // in other direction absent
        }
        throw forbiddenGraph;
    }
}

void ThicknessTwoCheckerMulti::checkProperty(const vector<adjacency_matrix_t> &matrices)
{
    auto m1 = matrices[1];

    auto edges = testPlanarity(m1);
    if (!edges.empty())
    {
        // printf("[");
        vector<forbidden_graph_t> forbiddenGraphs(3);
        for (auto e : edges)
        {
            // printf("(%d,%d),", e.first, e.second);
            forbiddenGraphs[1].push_back(make_pair(truth_value_true, e));
        }
        // printf("]\n");
        throw forbiddenGraphs;
    }

    auto m2 = matrices[2];
    edges = testPlanarity(m2);
    if (!edges.empty())
    {
        // printf("[");
        vector<forbidden_graph_t> forbiddenGraphs(3);
        for (auto e : edges)
        {
            // printf("(%d,%d),", e.first, e.second);
            forbiddenGraphs[2].push_back(make_pair(truth_value_true, make_pair(max(e.first, e.second), min(e.first, e.second)))); // from higher to lower vertex present
        }
        // printf("]\n");
        throw forbiddenGraphs;
    }
}
