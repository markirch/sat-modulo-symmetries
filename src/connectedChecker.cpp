#include "connectedChecker.hpp"

using std::make_pair;

// only works up to 32 vertices

void ConnectedChecker::checkProperty(const adjacency_matrix_t &matrix)
{
  int n = (int)matrix.size();
  // TODO: find connected components, ensure that at least one edge between them.

  vector<bool> reached(n, false);
  reached[0] = true;
  vector<int> expand = {0};

  while (!expand.empty())
  {
    int v = expand.back();
    expand.pop_back();
    for (int i = 0; i < n; i++)
    {
      if (matrix[v][i] == truth_value_true && !reached[i])
      {
        expand.push_back(i);
        reached[i] = true;
      }
    }
  }

  for (int i = 1; i < n; i++)
  {
    if (!reached[i])
    {
      // disconnected
      forbidden_graph_t forbidden_graph;
      // create reason
      for (int v1 = 0; v1 < n; v1++)
        for (int v2 = v1 + 1; v2 < n; v2++)
        {
          if (reached[v1] != reached[v2])
            forbidden_graph.push_back(make_pair(truth_value_false, make_pair(v1, v2)));
        }

      throw forbidden_graph;
    }
  }
}

bool checkConnected(const adjacency_matrix_t &matrix, const vector<bool> removedVertices, vector<int> &component, const int numRemoved)
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
  vector<bool> visited(n, false);
  vector<int> queue;
  visited[minVertex] = true;
  int numVisited = 1;
  queue.push_back(minVertex);
  while (!queue.empty())
  {
    int v = queue.back();
    queue.pop_back();
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

  if (numVisited < n - numRemoved)
  {
    for (int i = 0; i < n; i++)
      if (visited[i] && !removedVertices[i])
        component.push_back(i);
    return false;
  }

  return true;
}

bool testKConnected(const int k, const adjacency_matrix_t &matrix, vector<int> &component, vector<bool> &removedVertices)
{
  int n = (int)matrix.size();
  if (k == 3)
    for (int i1 = 0; i1 < n; i1++)
      for (int i2 = i1 + 1; i2 < n; i2++)
      {

        removedVertices = vector<bool>(n, false);
        removedVertices[i1] = true;
        removedVertices[i2] = true;
        if (!checkConnected(matrix, removedVertices, component, k - 1))
        {
          return false;
        }
      }
  else if (k == 2)
    for (int i1 = 0; i1 < n; i1++)
    {
      removedVertices = vector<bool>(n, false);
      removedVertices[i1] = true;
      if (!checkConnected(matrix, removedVertices, component, k - 1))
      {
        return false;
      }
    }
  else if (k == 1)
  {
    removedVertices = vector<bool>(n, false);
    if (!checkConnected(matrix, removedVertices, component, k - 1))
    {
      return false;
    }
  }
  else
  {
    printf("Error: k = %d not supported\n", k);
    EXIT_UNWANTED_STATE
  }
  return true;
}

void KConnectedChecker::checkProperty(const adjacency_matrix_t &matrix)
{
  vector<int> component;
  vector<bool> removedVertices;
  if (!testKConnected(k, matrix, component, removedVertices))
  {
    forbidden_graph_t forbiddenGraph;
    int n = (int)matrix.size();
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
    throw forbiddenGraph;
  }
}
