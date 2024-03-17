#include "forbiddenSubgraph.hpp"

HomomorphismParams defaultHomParams = {
    .timeout = std::make_shared<Timeout>(0s),
    .restarts_schedule = std::make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier),
    .no_supplementals = true};

InputGraph toGlasgowGraph(const adjacency_matrix_t &adjacencyMatrix, int nVertices)
{
    InputGraph G(nVertices, false, false);
    for (int i = 0; i < nVertices; i++)
        for (int j = i + 1; j < nVertices; j++)
            if (adjacencyMatrix[i][j] == truth_value_true)
                G.add_edge(i, j);
    return G;
}

InputGraph toLabeledGlasgowGraph(const adjacency_matrix_t &adjacencyMatrix, int nVertices)
{
    InputGraph G(nVertices, false, true); // have to mark that edge labels are present
    for (int i = 0; i < nVertices; i++)
        for (int j = i + 1; j < nVertices; j++)
        {
            if (adjacencyMatrix[i][j] == truth_value_true)
            {
                G.add_directed_edge(i, j, PRESENT_LABEL);
                G.add_directed_edge(j, i, PRESENT_LABEL);
            }
            else if (adjacencyMatrix[i][j] == truth_value_false)
            {
                G.add_directed_edge(i, j, ABSENT_LABEL);
                G.add_directed_edge(j, i, ABSENT_LABEL);
            }
            else
            {
                G.add_directed_edge(i, j, UNKNOWN_LABEL);
                G.add_directed_edge(j, i, UNKNOWN_LABEL);
            }
        }
    return G;
}

void HomomorphismResultToSubgraph(HomomorphismResult &res, forbidden_graph_t &forbidden_graph, const InputGraph &H)
{
    auto map_and_forbid_edge = [&res, &forbidden_graph](int u, int v, std::string_view )
    {
        if (u < v) // TODO directed version
        {
            forbidden_graph.push_back(std::make_pair(truth_value_true, std::make_pair(res.mapping[u], res.mapping[v])));
        }
    };
    H.for_each_edge(map_and_forbid_edge);
}

void InducedHomomorphismResultToSubgraph(HomomorphismResult &res, forbidden_graph_t &forbidden_graph, const InputGraph &H)
{
    int n = H.size();
    // printf("Size of mapping: %ld\n %d %d ,%d %d, %d %d \n", res.mapping.size(), 0, res.mapping[0], 1, res.mapping[1], 2, res.mapping[2]);

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            int pi = res.mapping[i]; // check whether correct direction
            int pj = res.mapping[j];
            if (H.edge_label(i, j) == PRESENT_LABEL)
                forbidden_graph.push_back(std::make_pair(truth_value_true, std::make_pair(pi, pj)));
            else // (H.edge_label(i, j) == ABSENT_LABEL)
                forbidden_graph.push_back(std::make_pair(truth_value_false, std::make_pair(pi, pj)));
            // else
            // {
            //     EXIT_UNWANTED_STATE
            // }
        }
}

// remove low degree vertices and check connected components
std::vector<std::vector<bool>> preprocess(const std::vector<std::vector<bool>> &G)
{
    int n = (int)G.size();
    if (n < 10)
        return G;
    // printf("G %ld\n", G.size());
    // for (int i = 0; i < n; i++)
    //     printf("\t%ld\n", G[i].size());
    for (int i = 0; i < n; i++)
    {

        int degree;
        degree = 0;
        for (int j = 0; j < n; j++)
        {
            if (G[i][j])
                degree++;
        }

        if (degree <= 2)
        {
            std::vector<std::vector<bool>> Gnew(G.begin(), G.begin() + n - 1);
            for (int r = 0; r < n - 1; r++)
            {
                Gnew[r].pop_back();
            }

            if (i < n - 1)
            {

                for (int j = 0; j < n - 1; j++)
                {
                    if (j == i)
                        continue;
                    Gnew[i][j] = Gnew[i][j] = G[n - 1][j];
                }
            }
            return preprocess(Gnew);
        }
    }

    return G;
}

void ForbiddenSubgraphCheckerGlasgow::checkProperty(const adjacency_matrix_t &matrix)
{
    // printAdjacencyMatrix(matrix, true);
    // int n = (int) matrix.size();
    // std::vector<std::vector<bool>> G1(n, std::vector<bool>(n, false));
    // for (int i = 0; i < n; i++)
    //     for (int j = i + 1; j < n; j++)
    //         if (matrix[i][j] == truth_value_true)
    //             G1[i][j] = G1[j][i] = true;

    // auto G2 = preprocess(G1);
    // if (G2.size() < 10)
    //     return;

    InputGraph G = toGlasgowGraph(matrix, (int)matrix.size());
    for (const InputGraph &H : forbiddenSubgraphsGlasgow)
    {
        int Gm = G.number_of_directed_edges();
        int Hm = H.number_of_directed_edges();
        if (Gm < Hm)
        {
            continue;
        }
        HomomorphismResult res = solve_homomorphism_problem(H, G, defaultHomParams);
        if (!res.mapping.empty())
        {
            forbidden_graph_t forbidden_graph;
            HomomorphismResultToSubgraph(res, forbidden_graph, H);
            throw forbidden_graph;
        }
    }

    InputGraph G2 = toLabeledGlasgowGraph(matrix, (int)matrix.size());
    for (const InputGraph &H : forbiddenInducedSubgraphsGlasgow)
    {
        // int Gm = G.number_of_directed_edges();
        // int Hm = H.number_of_directed_edges();
        // if (Gm < Hm)
        // {
        //     continue;
        // }
        HomomorphismResult res = solve_homomorphism_problem(H, G2, inducedHomParams);
        if (!res.mapping.empty())
        {
            forbidden_graph_t forbidden_graph;
            InducedHomomorphismResultToSubgraph(res, forbidden_graph, H);
            throw forbidden_graph;
        }
    }
}
