#include "connectedChecker.hpp"

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
