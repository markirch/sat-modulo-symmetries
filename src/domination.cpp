#include "domination.hpp"
#include <cassert>

bool checkConnected(const int sizeOfInterface, const int numInterfaces, const adjacency_matrix_t &matrix, const vector<bool> removedVertices, vector<int> &component)
{
    // check if the graph is connected if not component contains the vertices of one connected component
    int n = (int)matrix.size();
    int minVertex = -1; // smallest vertex in the graph which wasn't removed
    for (int i = 0; i < n; i++)
    {
        if (!removedVertices[i])
        {
            minVertex = i;
            break;
        }
    }
    assert(minVertex != -1);
    vector<bool> visited(n, false);
    vector<int> queue;
    visited[minVertex] = true;
    int numVisited = 1;
    queue.push_back(minVertex);
    while (!queue.empty())
    {
        int v = queue.back();
        queue.pop_back();

        if (v < numInterfaces * sizeOfInterface) // is a vertex in the interface
        {
            int interfaceIndex = v / sizeOfInterface;
            // mark all other vertices in the interface
            for (int j2 = 0; j2 < sizeOfInterface; j2++)
            {
                int vertex2 = interfaceIndex * sizeOfInterface + j2;
                if (!visited[vertex2] && !removedVertices[vertex2])
                {
                    visited[vertex2] = true;
                    numVisited++;
                    queue.push_back(vertex2);
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            if (matrix[v][i] == truth_value_true && !visited[i] && !removedVertices[i])
            {
                visited[i] = true;
                numVisited++;
                queue.push_back(i);
            }
        }
    }

    if (numVisited < n - (sizeOfInterface - 1))
    {
        for (int i = 0; i < n; i++)
            if (visited[i] && !removedVertices[i])
                component.push_back(i);
        return false;
    }

    return true;
}

bool specialConnectedness(const int sizeOfInterface, const int numInterfaces, const adjacency_matrix_t &matrix, vector<int> &component, vector<bool> &removedVertices)
{
    // try deleting k - 1 vertices and check if still connected assuming all vertices are additionally connected by a vertex.

    assert(sizeOfInterface <= 3 && sizeOfInterface >= 1);
    int n = (int)matrix.size();

    if (sizeOfInterface == 3)
        for (int i1 = 0; i1 < n; i1++)
            for (int i2 = i1 + 1; i2 < n; i2++)
            {

                removedVertices = vector<bool>(n, false);
                removedVertices[i1] = true;
                removedVertices[i2] = true;
                if (!checkConnected(sizeOfInterface, numInterfaces, matrix, removedVertices, component))
                {
                    return false;
                }
            }
    else if (sizeOfInterface == 2)
        for (int i1 = 0; i1 < n; i1++)
        {
            removedVertices = vector<bool>(n, false);
            removedVertices[i1] = true;
            if (!checkConnected(sizeOfInterface, numInterfaces, matrix, removedVertices, component))
            {
                return false;
            }
        }
    else if (sizeOfInterface == 1)
    {
        removedVertices = vector<bool>(n, false);
        if (!checkConnected(sizeOfInterface, numInterfaces, matrix, removedVertices, component))
        {
            return false;
        }
    }
    return true;
}

void QuasiKConnectedPropagator::checkProperty(const adjacency_matrix_t &matrix)
{
    vector<int> component;
    vector<bool> removedVertices;
    int n = (int)matrix.size();
    if (!specialConnectedness(interfaceSize, numInterfaces, matrix, component, removedVertices))
    {
        // the graph is not quasi k-connected
        // add the clause that at least one of the vertices in the component is not in the interface
        forbidden_graph_t forbiddenGraph;
        for (int i = 0; i < n; i++)
        {
            if (removedVertices[i])
                continue;

            for (int j = 0; j < n; j++)
            {
                if (removedVertices[j])
                    continue;

                // if exactly one of them is in the component
                if ((find(component.begin(), component.end(), i) != component.end()) != (find(component.begin(), component.end(), j) != component.end()))
                    forbiddenGraph.push_back({truth_value_false, {i, j}});
            }
        }

        // number of removed vertices
        int numRemoved = 0;
        for (int i = 0; i < n; i++)
            if (removedVertices[i])
                numRemoved++;

        // check that either all vertices of an interface are in the component or none (unless removed)
        for (int i = 0; i < numInterfaces; i++)
        {
            int numInComponent = 0;
            int numRemoved = 0;
            for (int j = 0; j < interfaceSize; j++)
            {
                if (removedVertices[i * interfaceSize + j])
                {
                    numRemoved++;
                    continue;
                }
                int vertex = i * interfaceSize + j;
                if (find(component.begin(), component.end(), vertex) != component.end())
                    numInComponent++;
            }
            if (numInComponent != 0 && numInComponent != interfaceSize - numRemoved)
            {
                printf("Removed vertices: ");
                for (int j = 0; j < n; j++)
                    if (removedVertices[j])
                        printf("%d ", j);
                printf("\n");
                printf("Component: ");
                for (int j = 0; j < (int)component.size(); j++)
                    printf("%d ", component[j]);
                printf("\n");
                EXIT_UNWANTED_STATE
            }
        }

        throw forbiddenGraph;
    }
};