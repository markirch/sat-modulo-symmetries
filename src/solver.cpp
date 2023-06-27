/**
 * Implementation of common parts amongst all solvers
 */
#include "solver.hpp"
#include <cmath>

// returns false if a clause was added
bool GraphSolver::cutoffFunction()
{
  if (config.assignmentCutoffPrerun && config.assignmentCutoffPrerun > stats.callsPropagator)
    return true;

  if (config.assignmentCutoffPrerunTime && config.assignmentCutoffPrerunTime > (clock() - stats.start) / CLOCKS_PER_SEC)
    return true;

  config.assignmentCutoffPrerunTime = 0;
  config.assignmentCutoffPrerun = 0;
  // frequency = 1; // check each time when creating cubes to add as much literals as possible

  adjacency_matrix_t matrix = getAdjacencyMatrix();

  if (config.assignmentScoring == COUNT_ASSIGNED)
  {
    int nAssigned = 0;
    for (int i = 0; i < vertices; i++)
      for (int j = i + 1; j < vertices; j++)
        if (matrix[i][j] != truth_value_unknown)
          nAssigned++;

    if (nAssigned <= config.assignmentCutoff)
      return true;
  }
  else
  {
    // score = sum log_2 ( 1 / p_e )
    // where p_e is the empirically recorded frequency of e
    // for balanced variables that are true half the time,
    // the contribution towards the score is 1
    // this is supposed to measure the search-space reduction
    // in the number of bits
    double assignmentScore = 0.0;
    for (int i = 0; i < vertices; i++)
    {
      for (int j = i + 1; j < vertices; j++)
      {
        double p_inv = stats.callsCheck;
        if (matrix[i][j] == truth_value_true)
        {
          if (edge_stats[i][j])
          {
            p_inv = p_inv / edge_stats[i][j];
          }
          assignmentScore += log2(p_inv);
        }
        else if (matrix[i][j] == truth_value_false)
        {
          if (edge_stats[i][j] != stats.callsCheck)
          {
            p_inv = p_inv / (stats.callsCheck - edge_stats[i][j]);
          }
          assignmentScore += log2(p_inv);
        }
      }
    }

    if (assignmentScore <= config.assignmentCutoff)
    {
      return true;
    }
  }

  if (!checkPartiallyDefined(true)) // check all before creating cube
    return false;

  for (int i = 0; i < vertices; i++)
    for (int j = i + 1; j < vertices; j++)
    {
      if (matrix[i][j] == truth_value_true)
        printf(" %d", edges[i][j]);
      if (matrix[i][j] == truth_value_false)
        printf(" -%d", edges[i][j]);
    }
  printf("\n");
  vector<lit_t> clause;
  for (int i = 0; i < vertices; i++)
    for (int j = i + 1; j < vertices; j++)
    {
      if (matrix[i][j] == truth_value_true)
        clause.push_back(-edges[i][j]);
      if (matrix[i][j] == truth_value_false)
        clause.push_back(edges[i][j]);
    }

  // printf("Size %ld\n", clause.size());
  addClause(clause, false);
  return false;
}

bool GraphSolver::propagate()
{
  stats.callsPropagator++;

  if (config.printPartiallyDefined && rand() % config.printPartiallyDefined == 0)
  {
    auto matrix = getAdjacencyMatrix();
    printf("p\t");
    printPartiallyDefinedAdjacencyMatrix(matrix);
  }

  auto start = clock();
  bool res = true;
  auto matrix = getAdjacencyMatrix();
  res = checkPartiallyDefined(false);

  if (res && config.assignmentCutoff)
  {
    res = cutoffFunction();
  }

  // if (checkEmbeddabilityKS && rand() % checkEmbeddabilityKS == 0)
  // {
  //   auto m = getAdjacencyMatrix();
  //   res = testEmebeddabilityKS(m);
  //   if (!res)
  //   {
  //     clause_t clause;
  //     for (int i = 0; i < m.size(); i++)
  //     {
  //       for (int j = i + 1; j < m.size(); j++)
  //       {
  //         if (m[i][j] == truth_value_true)
  //           clause.push_back(-edges[i][j]);
  //       }
  //     }
  //     addClause(clause, false);
  //   }
  // }

  // if (rand() % 10000 == 0)
  //  res = res && extractForbiddenSubgraphs(getAdjacencyMatrix());

  /*if (res && hypercoloring) {
    adjacency_matrix_t matrix = getAdjacencyMatrix();
    bool sufficiently_populated = true;
    for (int e = b_vertices[0]; e < vertices; e++) {
      int e_size = 0;
      for (int v = 0; v < b_vertices[0]; v++) {
        if (matrix[v][e] == truth_value_true) {
          e_size++;
          if (e_size == 2) {
            break;
          }
        }
      }
      if (e_size < 2) {
        sufficiently_populated = false;
        break;
      }
    }
    if (sufficiently_populated) {
      coloring_t coloring(b_vertices[0]);
      if (getHyperColoring(b_vertices, matrix, coloring)) {
        res = false;
        addClause(getHyperColoringClause(coloring, b_vertices, matrix, edges), false);
        stats.hyperclauses++;
      }
    }
  }*/

  stats.timePropagator += clock() - start;
  return res;
}

bool GraphSolver::checkPartiallyDefined(bool isFullyDefined)
{
  adjacency_matrix_t matrix = getAdjacencyMatrix();
  try
  {
    for (auto checker : this->partiallyDefinedGraphCheckers)
      checker->check(matrix, isFullyDefined);
  }
  catch (const forbidden_graph_t forbiddenGraph)
  {
    // trnsform forbidden subgraph into a clause, which blocks this graph
    clause_t clause;
    for (auto signedEdge : forbiddenGraph)
    {
      auto edge = signedEdge.second;
      if (signedEdge.first == truth_value_true)
        clause.push_back(-edges[edge.first][edge.second]);
      else // assum that not truth_value_unknown
        clause.push_back(edges[edge.first][edge.second]);
    }
    addClause(clause, false); // TODO eventually vector with checkers which are redundant and which are not; only makes sense when supported by Cadical
    return false;
  }

  try
  {
    vector<truth_value_t> curAssignment = getCurrentAssignemnt();
    for (auto checker : this->complexPartiallyDefinedGraphCheckers)
      checker->check(matrix, curAssignment, isFullyDefined);
  }
  catch (const vector<clause_t> clauses)
  {
    for (clause_t clause : clauses)
      addClause(clause, false);
    return false;
  }

  return true;
}

// returns false if graph does not satisfy a property
bool GraphSolver::checkFullyDefinedGraph(const adjacency_matrix_t &matrix, const vector<int> &model)
{
  try
  {
    for (auto checker : this->fullyDefinedGraphCheckers)
      checker->check(matrix);
  }
  catch (const forbidden_graph_t forbiddenGraph)
  {
    // trnsform forbidden subgraph into a clause, which blocks this graph
    clause_t clause;
    for (auto signedEdge : forbiddenGraph)
    {
      auto edge = signedEdge.second;
      if (signedEdge.first == truth_value_true)
        clause.push_back(-edges[edge.first][edge.second]);
      else // assum that not truth_value_unknown
        clause.push_back(edges[edge.first][edge.second]);
    }
    addClause(clause, false); // TODO eventually vector with checkers which are redundant and which are not; only makes sense when supported by Cadical
    return false;
  }

  try
  {
    vector<truth_value_t> curAssignment = getCurrentAssignemnt();
    for (auto checker : this->complexFullyDefinedGraphCheckers)
      checker->check(matrix, model, nextFreeVariable);
  }
  catch (const vector<clause_t> clauses)
  {
    for (clause_t clause : clauses)
      addClause(clause, false);
    return false;
  }

  return true; // all checks passed
}

bool GraphSolver::check()
{
  if (!checkPartiallyDefined(true))
    return false;

  // ----------------- all checks done for partially defined graphs -------------------------------------------
  // printf("Start fully defined\n");
  adjacency_matrix_t matrix = getAdjacencyMatrix();
  stats.callsCheck++;
  if (config.printFullyDefinedGraphs)
  {
    printf("Check\n");
    fflush(stdout);
    if (config.hypermode)
    {
      printHypergraph(matrix, config.b_vertices);
    }
    else
    {
      printAdjacencyMatrix(matrix, config.printFullMatrix);
    }
  }
  recordGraphStats(matrix);
  // fully defined and minimal graph
  if (config.printIntermediateStatistic)
  {
    printf("Fully defined: %lld\n", stats.callsCheck);
    printStatistics();
  }

  clock_t start = clock();
  bool r = checkFullyDefinedGraph(matrix, *model);
  stats.timeCheckFullGraphs += clock() - start;
  if (!r)
    return false;

  nModels++;
  if (!config.hideGraphs && !config.quiet)
  {

    printf("Solution %d\n", nModels);
    if (config.hypermode)
    {
      printHypergraph(matrix, config.b_vertices);
    }
    else
    {
      printAdjacencyMatrix(matrix, config.printFullMatrix);
    }
  }

  if (config.allModels)
  {
    // exclude current graph
    adjacency_matrix_t &m = matrix;
    vector<lit_t> clause;

    if (config.hyperedgeColoring)
    {
      // forbid not just the specific hypergraph, but also its intersection graph
      adjacency_matrix_t im = getIntersectionMatrix(matrix, config.b_vertices);
      for (int i = 0; i < config.b_vertices[1]; i++)
        for (int j = i + 1; j < config.b_vertices[1]; j++)
          if (im[i][j] == truth_value_true)
            clause.push_back(-config.edges_intersection_graph[i][j]);
          else if (im[i][j] == truth_value_false)
            clause.push_back(config.edges_intersection_graph[i][j]);
          else
            EXIT_UNWANTED_STATE
      addClause(clause, false);
      clause.resize(0);
    }

#ifndef DIRECTED
    for (int i = 0; i < vertices; i++)
      for (int j = i + 1; j < vertices; j++)
        if (m[i][j] == truth_value_true)
          clause.push_back(-edges[i][j]);
        else if (m[i][j] == truth_value_false)
          clause.push_back(edges[i][j]);
        else
          EXIT_UNWANTED_STATE
    addClause(clause, false);
    return false;
#else
    for (int i = 0; i < vertices; i++)
      for (int j = 0; j < vertices; j++)
        if (i == j)
        {
        }
        else if (m[i][j] == truth_value_true)
          clause.push_back(-edges[i][j]);
        else if (m[i][j] == truth_value_false)
          clause.push_back(edges[i][j]);
        else
          EXIT_UNWANTED_STATE
    addClause(clause, false);
    return false;
#endif
  }
  return true;
}

void GraphSolver::printStatistics()
{
  if (!config.printStats)
    return;
  printf("Time in propagator: %f\n", ((double)stats.timePropagator) / CLOCKS_PER_SEC);
  printf("Time in check full graphs: %f\n", ((double)stats.timeCheckFullGraphs) / CLOCKS_PER_SEC);
  printf("Calls of check: %lld\n", stats.callsCheck);
  printf("Calls propagator: %lld\n", stats.callsPropagator);
  if (config.allModels)
    printf("Number of models: %d\n", nModels);
  // printEdgeStats();

  vector<GraphChecker *> allCheckers;
  allCheckers.insert(allCheckers.end(), partiallyDefinedGraphCheckers.begin(), partiallyDefinedGraphCheckers.end());
  allCheckers.insert(allCheckers.end(), complexPartiallyDefinedGraphCheckers.begin(), complexPartiallyDefinedGraphCheckers.end());
  allCheckers.insert(allCheckers.end(), fullyDefinedGraphCheckers.begin(), fullyDefinedGraphCheckers.end());
  allCheckers.insert(allCheckers.end(), complexFullyDefinedGraphCheckers.begin(), complexFullyDefinedGraphCheckers.end());
  for (auto checker : allCheckers)
    checker->printStats();
}

void GraphSolver::solve()
{
  // solve
  if (!config.quiet)
  {
    printf("Start solving\n");
    fflush(stdout);
  }

  // get a solve handle
  int cubeCounter = 0;
  if (!config.cubeFile.empty())
  {
    ifstream is(config.cubeFile);
    string line;
    while (getline(is, line))
    {
      cubeCounter++;
      if (config.rangeCubes.first != 0)
      {
        if (cubeCounter < config.rangeCubes.first)
          continue;
        if ((cubeCounter > config.rangeCubes.second))
          break;
      }

      initEdgeMemory();
      if (config.non010colorable)
        initTriangleMemory();
      printf("Solve cube %d\n", cubeCounter);
      vector<lit_t> assumptions; // variable names from the initial encoding

      printf("%s\n", line.c_str());

      istringstream iss(line);

      string lit;
      while (std::getline(iss, lit, ' '))
      {
        if (lit.empty())
          continue;
        assumptions.push_back(stoi(lit));
      }

      clock_t start = clock();

      bool solvedSuccessfully = true;
      if (config.timeout)
      {
        if (!solve(assumptions, config.timeout))
          printf("Timeout reached\n");
      }
      else
      {
        solve(assumptions);
      }
      if (solvedSuccessfully)
        printf("Time for cube %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
      printStatistics();
    }
    printf("All cubes solved\n");
  }
  else
  {
    initEdgeMemory();
    if (config.non010colorable)
      initTriangleMemory();
    if (config.timeout)
    {
      if (!solve(vector<int>(), config.timeout))
        printf("Timeout reached\n");
    }
    else
    {
      solve(vector<int>());
    }

    if (!config.quiet)
    {
      printf("Searched finished\n");
      printStatistics();
    }
  }
}

void GraphSolver::initEdgeMemory()
{
  edge_stats.resize(vertices);
  for (int v = 0; v < vertices; v++)
    edge_stats[v].resize(vertices, 0);
}

void GraphSolver::initTriangleMemory()
{
  triangle_stats.resize(vertices);
  for (int v = 0; v < vertices; v++)
  {
    triangle_stats[v].resize(vertices);
    for (int w = 0; w < vertices; w++)
    {
      triangle_stats[v][w].resize(vertices, 0);
    }
  }
}

void GraphSolver::printEdgeStats()
{
  printf("Edge occurrence statistics:\n");
  for (int u = 0; u < vertices; u++)
  {
    for (int v = 0; v < vertices; v++)
      printf("%4u ", edge_stats[u][v]);
    printf("\n");
  }
}

/*
static inline void greedyConstraint(clingo_propagate_control_t *ctl, int *vertexColoring, propagator_t *data)
{
  nColorings++;
  // printf("START\n");
  int colors = 4;
  int n = data->nVertices;
  int vertexOrdering[n];

  int pos = 0;
  for (int c = 0; c < colors; c++)
    for (int i = 0; i < n; i++)
      if (vertexColoring[i] == c)
        vertexOrdering[pos++] = i;

  assert(pos == n);

  // coloring[i][0] denotes if vertex vertexOrdering[i] has color 0
  clingo_literal_t coloring[n][colors];  // coloring of each vertex
  clingo_literal_t available[n][colors]; // check if color is available, i.e., no smaller vertex in ordering has the color
  clingo_literal_t isColored[n];         // check if vertex is colored
  clingo_literal_t adjacentAndColorC[n][n][colors];

  for (int i = 0; i < n; i++)
  {
    if (!clingo_propagate_control_add_literal(ctl, &isColored[i]))
      exit(EXIT_FAILURE);

    for (int c = 0; c < colors; c++)
    {
      if (!clingo_propagate_control_add_literal(ctl, &coloring[i][c]))
        exit(EXIT_FAILURE);
      if (!clingo_propagate_control_add_literal(ctl, &available[i][c]))
        exit(EXIT_FAILURE);

      for (int j = i + 1; j < n; j++)
      {
        if (!clingo_propagate_control_add_literal(ctl, &adjacentAndColorC[i][j][c]))
          exit(EXIT_FAILURE);

        // printf("%d\n", adjacentAndColorC[i][j][c]);
      }
    }
  }

  // --------------add clauses-----------

  bool res;

  // one color implies isColored
  for (int i = 0; i < n; i++)
  {
    for (int c = 0; c < colors; c++)
    {
      clingo_literal_t clause[] = {-coloring[i][c], isColored[i]};
      addClauseToBuff(clause, 2);
    }
  }

  // at most on color
  for (int i = 0; i < n; i++)
  {
    for (int c1 = 0; c1 < colors; c1++)
    {
      for (int c2 = c1 + 1; c2 < colors; c2++)
      {
        clingo_literal_t clause[] = {-coloring[i][c1], -coloring[i][c2]};
        addClauseToBuff(clause, 2);
      }
    }
  }

  // at least one node is not colored
  clingo_literal_t clause[n];
  for (int i = 0; i < n; i++)
    clause[i] = -isColored[i];
  addClauseToBuff(clause, n);

  // smallest available means that it should get this color
  for (int i = 0; i < n; i++)
  {
    for (int c1 = 0; c1 < colors; c1++)
    {
      clingo_literal_t clause[colors + 2];
      for (int c2 = 0; c2 < c1; c2++) // smaller colors
        clause[c2] = available[i][c2];

      clause[c1] = -available[i][c1];
      clause[c1 + 1] = coloring[i][c1];

      addClauseToBuff(clause, c1 + 2);
    }
  }

  // truth_value_true and color c
  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      for (int c = 0; c < colors; c++)
      {

        clingo_literal_t clause1[] = {coloring[i][c], -adjacentAndColorC[i][j][c]}; // not truth_value_true implies not truth_value_true and color c
        addClauseToBuff(clause1, 2);

        int v1 = MIN(vertexOrdering[i], vertexOrdering[j]);
        int v2 = MAX(vertexOrdering[i], vertexOrdering[j]);
        clingo_literal_t clause2[2] = {data->edgeLits[v1][v2], -adjacentAndColorC[i][j][c]}; // not truth_value_true implies not truth_value_true and color c
        addClauseToBuff(clause2, 2);

        clingo_literal_t clause3[] = {-data->edgeLits[v1][v2], -coloring[i][c], adjacentAndColorC[i][j][c]};
        addClauseToBuff(clause3, 3);
      }
    }
  }

  // if no smaller vertex is truth_value_true and has color c, than color c is available.
  for (int i = 0; i < n; i++)
  {
    clingo_literal_t clause[i];
    for (int c = 0; c < colors; c++)
    {
      for (int j = 0; j < i; j++)
        clause[j] = adjacentAndColorC[j][i][c];
      clause[i] = available[i][c];
      addClauseToBuff(clause, i + 1);
    }
  }

  // if truth_value_true and color c than not available
  for (int i = 0; i < n; i++)
  {
    for (int c = 0; c < colors; c++)
    {
      for (int j = 0; j < i; j++)
      {
        // TODO maybe use edge literals and color literals directly.
        clingo_literal_t clause[] = {-adjacentAndColorC[j][i][c], -available[i][c]};
        addClauseToBuff(clause, 2);
      }
    }
  }
  // fflush(stdout);
} */

void GraphSolver::recordGraphStats(const adjacency_matrix_t &matrix)
{
  for (int u = 0; u < vertices; u++)
  {
    for (int v = 0; v < vertices; v++)
    {
      if (matrix[u][v])
      {
        edge_stats[u][v]++;
        if (config.non010colorable)
        {
          for (int w = 0; w < vertices; w++)
          {
            if (matrix[u][w] && matrix[v][w])
            {
              triangle_stats[u][v][w]++;
            }
          }
        }
      }
    }
  }
}
