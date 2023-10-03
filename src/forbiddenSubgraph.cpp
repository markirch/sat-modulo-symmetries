#include "forbiddenSubgraph.hpp"

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

#ifdef GLASGOW

void ForbiddenSubgraphCheckerGlasgow::checkProperty(const adjacency_matrix_t &matrix)
{
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
            auto map_and_forbid_edge = [&res, &forbidden_graph](int u, int v, std::string_view s)
            {
                if (u < v) // TODO directed version
                {
                    forbidden_graph.push_back(std::make_pair(truth_value_true, std::make_pair(res.mapping[u], res.mapping[v])));
                }
            };
            H.for_each_edge(map_and_forbid_edge);
            throw forbidden_graph;
        }
    }
}

#endif