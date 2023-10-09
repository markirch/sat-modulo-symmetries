#include "minimalityCheck.hpp"
#include "useful.h"

using std::pair;
using std::make_pair;

/**
 * @brief checks whether swapping connected components of two colors can decrease the adjacency matrix, unknown edges are assumed to be present for the connected componenets
 *
 * @param matrix the adjacency matrix
 */
bool checkSwaps(adjacency_matrix_t &matrix, vector<pair<int, int>> intervallsColoring, FILE *symBreakClausesFile)
{
    int N = intervallsColoring.size(); // last partition is only the single vertex
    int vertices = (int)matrix.size();
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++)
        {
            // get connected components from vertices form color i and j

            vector<int> componenet(vertices, 0); // 0 if not assigned to a component yet

            // get unmarked vertex
            int unmarkedVertex = intervallsColoring[i].first;
            int nComponents = 0;

            do
            {
                nComponents++;

                vector<int> expand;
                componenet[unmarkedVertex] = nComponents;
                expand.push_back(unmarkedVertex);

                while (!expand.empty())
                {
                    int v = expand[expand.size() - 1];
                    expand.pop_back();

                    for (int u = 0; u < vertices; u++)
                    {
                        if (u == v || componenet[u] == nComponents || matrix[v][u] == truth_value_false)
                            continue;

                        if ((u >= intervallsColoring[i].first && u < intervallsColoring[i].second) ||
                            (u >= intervallsColoring[j].first && u < intervallsColoring[j].second))
                        {
                            componenet[u] = nComponents;
                            expand.push_back(u);
                        }
                    }
                }

                // check number of vertices of each color in the connected componenet
                int ni = 0; // vertices with color i
                int nj = 0; // vertices with color j

                for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                    if (componenet[v] == nComponents)
                        ni++;
                for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                    if (componenet[v] == nComponents)
                        nj++;
                // printf("Size of component %d: %d %d\n", c, ni, nj);
                /* if (nj > ni)
                {
                    printf("Could be discarded %d %d; %d %d !!!\n", ni, nj, i, j);
                    return false;
                } */

                if (ni == nj && ni != intervallsColoring[i].second - intervallsColoring[i].first) // check connected components if not all vertices
                {
                    partition_t partition(vertices, false);
                    for (int i = 0; i < N; i++)
                    {
                        partition[intervallsColoring[i].first] = true;
                        partition[intervallsColoring[i].first + 1] = true; // first vertex in own parition.
                    }
                    partition[vertices - 1] = true; // last vertex in own parition

                    // check swapping this two blocks
                    vertex_ordering_t vertex_ordering(vertices);
                    for (int i = 0; i < vertices; i++)
                        vertex_ordering[i] = i;

                    int pos1 = intervallsColoring[i].first;
                    int pos2 = intervallsColoring[j].first;

                    for (int k = 0; k < ni; k++)
                    {
                        // find the next two vertices to swap (in each color seperately)
                        while (componenet[pos1] != nComponents)
                            pos1++;

                        while (componenet[pos2] != nComponents)
                            pos2++;

                        std::swap(vertex_ordering[pos1], vertex_ordering[pos2]);

                        pos1++;
                        pos2++;
                    }
                    try
                    {
                        checkMinimality(matrix, vertex_ordering, {.initial_partition = partition, .cutoff = 0});
                    }
                    catch (std::vector<signed_edge_t> e)
                    {
                        // printf("Not minimal by connected components\n");

                        if (symBreakClausesFile)
                        {
                            fprintf(symBreakClausesFile, "c;%ld;", e.size());

                            // print connected component
                            for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                                if (componenet[v] == nComponents)
                                    fprintf(symBreakClausesFile, "%d ", v);

                            for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                                if (componenet[v] == nComponents)
                                    fprintf(symBreakClausesFile, "%d ", v);
                            fprintf(symBreakClausesFile, ";");

                            // print the two parititions used
                            fprintf(symBreakClausesFile, "%d ; %d;", i, j);
                        }

                        // add information that connected component
                        for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                        {
                            if (componenet[v] != nComponents)
                                continue;
                            for (int u = intervallsColoring[j].first; u < intervallsColoring[j].second; u++)
                            {
                                if (componenet[u] == nComponents)
                                    continue;
                                e.push_back(make_pair(truth_value_true, make_pair(v, u)));
                            }
                        }

                        for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                        {
                            if (componenet[v] == nComponents)
                                continue;
                            for (int u = intervallsColoring[j].first; u < intervallsColoring[j].second; u++)
                            {
                                if (componenet[u] != nComponents)
                                    continue;
                                e.push_back(make_pair(truth_value_true, make_pair(v, u)));
                            }
                        }

                        throw e;
                    }
                }

                // search for unmarked vertex
                unmarkedVertex = -1;
                for (int v = intervallsColoring[i].first; v < intervallsColoring[i].second; v++)
                    if (componenet[v] == 0)
                        unmarkedVertex = v;
                for (int v = intervallsColoring[j].first; v < intervallsColoring[j].second; v++)
                    if (componenet[v] == 0)
                        unmarkedVertex = v;
            } while (unmarkedVertex != -1);
        }

    return true;
}
