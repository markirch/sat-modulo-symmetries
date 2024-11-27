#include "useful.h"
#include "cadicalSMS.hpp"
#include "cadical.hpp"
#include <algorithm>

void CadicalSolver::init(SolverConfig config, cnf_t &cnf)
{
    // The root-level of the trail is always there
    current_trail.push_back(std::vector<int>());

    // only_propagating = false;
    solver = new CaDiCaL::Solver();
    if (!config.proofFile.empty())
        if (!solver->trace_proof(config.proofFile.c_str()))
            EXIT_UNWANTED_STATE
    // if (!solver->set("check", 1))
    //     EXIT_UNWANTED_STATE
    // if (!solver->set("checkproof", 1))
    //     EXIT_UNWANTED_STATE

    // TODO currently inprocessing not supported by cadical with theory propagation
    // if (!solver->set("compact", 0))
    //     EXIT_UNWANTED_STATE
    // solver->configure("plain");
    if (!solver->configure("unsat"))
        EXIT_UNWANTED_STATE
    // if (!solver->set("chrono", 0))
    //    EXIT_UNWANTED_STATE

    // solver->set("lucky", 0);
    if (config.turnOffInprocessing)
        solver->set("inprocessing", 0);
    // solver->set("log", 1);
    // solver->set("debug", 1);

    // solver->set("probeint", 1);

    if (config.cadicalConfig != "") {
      std::istringstream iss(config.cadicalConfig);
      string param;
      while (iss >> param) {
        if (!solver->set_long_option(("--" + param).c_str())) {
          std::cerr << "invalid Cadical option '" << param << "'" << std::endl;
          EXIT_UNWANTED_STATE
        }
      }
    }

    // register propagator first
    solver->connect_external_propagator(this);

    for (auto v : config.observedVars)
    {
        solver->add_observed_var(v);
    }

    int highestObservedVariables = *max_element(config.observedVars.begin(), config.observedVars.end());
    currentAssignment = vector<truth_value_t>(highestObservedVariables + 1, truth_value_unknown); // must be created before adding clauses because already notification can happen
    isFixed = vector<bool>(highestObservedVariables + 1, false);

    literal2clausePos = vector<vector<int>>(highestObservedVariables + 1);
    literal2clauseNeg = vector<vector<int>>(highestObservedVariables + 1);

    // add clauses to solver
    for (auto clause : cnf)
    {
        if (clause.size() == 0)
            EXIT_UNWANTED_STATE

        for (auto lit : clause)
        {
            if (lit == 0)
                EXIT_UNWANTED_STATE
            solver->add(lit);
        }
        solver->add(0);

        if (config.addedClauses)
        {
            for (auto lit : clause)
                fprintf(config.addedClauses, "%d ", lit);
            fprintf(config.addedClauses, "0\n");
        }
    }

    initEdgeMemory();
}

// add formula and register propagator
CadicalSolver::CadicalSolver(SolverConfig config, cnf_t &cnf) : GraphSolver(config)
{
    init(config, cnf);
}

CadicalSolver::CadicalSolver(SolverConfig config) : GraphSolver(config)
{
    cnf_t cnf; // empty cnf
    init(config, cnf);
}

int* CadicalSolver::getNextGraph(vector<int> assumptions)
{
    // mandatory configuration for this API to make sense
    config.checkSolutionInProp = false;
    config.allModels = false;
    config.hideGraphs = true;

    while (true) {

      for (auto lit : assumptions)
          solver->assume(lit);
      if (config.addedClauses) {
          if (assumptions.size() != 0)
              printf("Warning: adding assumptions to clause file; only supported for a single assumption");
          for (auto l : assumptions)
              fprintf(config.addedClauses, "%d 0\n", l); // unit clauses
      }

      int res = solver->solve();
      if (res == 20) // no more solutions
          return nullptr;
      if (res != 10) { EXIT_UNWANTED_STATE } // TODO just to be sure for know

      adjacency_matrix_t matrix = getAdjacencyMatrix();
      vector<lit_t> clause; // to block the current solution

      if (check_solution()) { // true if no clause was added
          last_graph.clear();
          int m = 0;
          last_graph.push_back(m);
          for (int i = 0; i < vertices; i++) {
            for (int j = i+1; j < vertices; j++) {
              if (matrix[i][j] == truth_value_true) {
                m++;
                last_graph.push_back(i);
                last_graph.push_back(j);
                clause.push_back(-edges[i][j]);
              } else if (matrix[i][j] == truth_value_false) {
                clause.push_back(edges[i][j]);
              } else {
                EXIT_UNWANTED_STATE
              }
            }
          }
          last_graph[0] = m;
          addClause(clause, false);
          return &(last_graph[0]);
      }

    }
}

bool CadicalSolver::solve(vector<int> assumptions)
{
    do
    {
        for (auto lit : assumptions)
            solver->assume(lit);
        if (config.addedClauses)
        {
            if (assumptions.size() != 0)
                printf("Warning: adding assumptions to clause file; only supported for a single assumption");
            for (auto l : assumptions)
                fprintf(config.addedClauses, "%d 0\n", l); // unit clauses
        }

        int res = solver->solve();
        if (res == 20) // not
            return false;
        if (res != 10)
        {
            return false; // return false if wasn't solved completely
        }

        if (check_solution()) // true if no clause was added
            return true;
    } while (true);
    return true;
}

void CadicalSolver::printFullModel()
{
    printf("Model: ");
    for (int i = 1; i < nextFreeVariable; i++)
    {
        printf("%d ", solver->val(i));
    }
    printf("\n");
}

bool CadicalSolver::solve(vector<int> assumptions, int timeout)
{
    TimeoutTerminator t = TimeoutTerminator(timeout);
    solver->connect_terminator(&t);
    bool res = solve(assumptions);
    solver->disconnect_terminator();
    return res;
}

// PySMS API
extern "C" {
  int* next_solution(void* sms_solver) {
    return ((CadicalSolver*) sms_solver)->getNextGraph(vector<int>{});
  }
  void* create_solver(int vertices) {
    return new(std::nothrow) CadicalSolver(SolverConfig(vertices));
  }
  void destroy_solver(void* sms_solver) {
    delete (CadicalSolver*) sms_solver;
  }
  void add_literal(void* sms_solver, int lit) {
    ((CadicalSolver*)sms_solver)->solver->add(lit);
  }
}

// // TODO select parts which should be used
//     if (false)
//     {
//         int solutionCounter = 0;
//         vector<int> vertexOrderingColoringPrevious;
//         vector<int> color0prev;
//         vector<int> color1prev;
//         while (true) // iterate over all models
//         {
//             for (auto lit : assumptions)
//                 solver->assume(lit);
//             int res = solver->solve();
//             if (res == 20)
//             {
//                 printf("All models found\n");
//                 break;
//             }
//             else if (res == 10)
//             {
//                 stats.callsCheck++;
//                 clock_t start = clock();
//                 // satisfiable
//                 adjacency_matrix_t matrix = vector<vector<truth_value_t>>(vertices, vector<truth_value_t>(vertices, truth_value_unknown));
//                 for (int i = 0; i < vertices; i++)
//                 {
//                     for (int j = i + 1; j < vertices; j++)
//                     {
//                         matrix[i][j] = matrix[j][i] = solver->val(edges[i][j]) > 0 ? truth_value_true : truth_value_false;
//                     }
//                 }

//                 // printf("Model:");
//                 // for (int i = 7889; i < min(highestVariable, 7938 + 1); i++)
//                 //     printf("%d ", solver->val(i));
//                 // printf("\n");

//                 if (0)
//                 {

//                     // printAdjacencyMatrix(matrix);

//                     // if (vertexOrderingColoringPrevious.size() != 0)
//                     // {

//                     //     vector<int> testColoring = vector<int>(vertices, 0);
//                     //     for (int i = 0; i < vertices; i++)
//                     //     {
//                     //         int a0 = true; // color 0 available
//                     //         int a1 = true; // color 1 available
//                     //         for (int j = 0; j < i; j++)
//                     //             if (matrix[vertexOrderingColoringPrevious[i]][vertexOrderingColoringPrevious[j]] == truth_value_true && testColoring[j] == 0)
//                     //             {
//                     //                 a0 = false;
//                     //                 printf("Reason0 for %d: %d\n", vertexOrderingColoringPrevious[i], vertexOrderingColoringPrevious[j]);
//                     //             }
//                     //         for (int j = 0; j < i; j++)
//                     //             for (int k = j + 1; k < i; k++)
//                     //                 if (matrix[vertexOrderingColoringPrevious[i]][vertexOrderingColoringPrevious[j]] == truth_value_true && matrix[vertexOrderingColoringPrevious[i]][vertexOrderingColoringPrevious[k]] == truth_value_true && matrix[vertexOrderingColoringPrevious[j]][vertexOrderingColoringPrevious[k]] == truth_value_true && testColoring[j] == 1 && testColoring[k] == 1)
//                     //                 {
//                     //                     a1 = false;
//                     //                     printf("Reason1 for %d: %d %d\n", vertexOrderingColoringPrevious[i], vertexOrderingColoringPrevious[j], vertexOrderingColoringPrevious[k]);
//                     //                 }
//                     //         if (a0)
//                     //             testColoring[i] = 0;
//                     //         else if (a1)
//                     //             testColoring[i] = 1;
//                     //         else
//                     //             testColoring[i] = -1;
//                     //     }

//                     //     for (int i = 0; i < vertices; i++)
//                     //     {
//                     //         int colSAt = -1;
//                     //         if (solver->val(color0prev[i]) > 0)
//                     //             colSAt = 0;
//                     //         if (solver->val(color1prev[i]) > 0)
//                     //             colSAt = 1;
//                     //         printf("(%d,%d, %d), ", vertexOrderingColoringPrevious[i], colSAt, testColoring[i]);
//                     //     }
//                     //     printf("\n");
//                     // }

//                     // printf("ASDF\n");
//                     coloring_t coloring(vertices);
//                     if (get010Coloring(vertices, matrix, coloring))
//                     {
//                         vector<int> vertexOrderingColoring;
//                         for (int i = 0; i < vertices; i++)
//                             if (coloring[i] == 0)
//                                 vertexOrderingColoring.push_back(i);

//                         for (int i = 0; i < vertices; i++)
//                             if (coloring[i] == 1)
//                                 vertexOrderingColoring.push_back(i);
//                         vertexOrderingColoringPrevious = vertexOrderingColoring;

//                         // printf("Ordering:");
//                         // for (auto v : vertexOrderingColoring)
//                         //     printf("%d ", v);
//                         // printf("\n");

//                         vector<int> color0;
//                         vector<int> color1;
//                         vector<int> uncolored;

//                         for (int i = 0; i < vertices; i++)
//                         {
//                             color0.push_back(++highestVariable);
//                             color1.push_back(++highestVariable);
//                             uncolored.push_back(++highestVariable);
//                         }

//                         color0prev = color0;
//                         color1prev = color1;

//                         // printf("Color0:");
//                         // for (auto v : color0)
//                         //     printf("%d ", v);
//                         // printf("\n");

//                         // printf("Color1:");
//                         // for (auto v : color1)
//                         //     printf("%d ", v);
//                         // printf("\n");

//                         for (int i = 0; i < vertices; i++)
//                         {
//                             int available0 = color0[i]; // if available then give color 0
//                             int available1 = ++highestVariable;

//                             vector<int> adjacentAndColor0;

//                             for (int j = 0; j < i; j++)
//                             {
//                                 int var = ++highestVariable;
//                                 solver->add(edges[vertexOrderingColoring[i]][vertexOrderingColoring[j]]);
//                                 solver->add(-var);
//                                 solver->add(0);

//                                 solver->add(color0[j]);
//                                 solver->add(-var);
//                                 solver->add(0);

//                                 solver->add(-edges[vertexOrderingColoring[i]][vertexOrderingColoring[j]]);
//                                 solver->add(-color0[j]);
//                                 solver->add(var);
//                                 solver->add(0);

//                                 adjacentAndColor0.push_back(var);
//                             }

//                             // \lnot available0 iff or( adjacentAndColor0)
//                             for (auto l : adjacentAndColor0)
//                             {
//                                 solver->add(-l);
//                                 solver->add(-available0);
//                                 solver->add(0);
//                             }
//                             for (auto l : adjacentAndColor0)
//                                 solver->add(l);
//                             solver->add(available0);
//                             solver->add(0);

//                             vector<int> adjacentTriangleAndColor1;
//                             for (int j = 0; j < i; j++)
//                             {
//                                 for (int k = j + 1; k < i; k++)
//                                 {
//                                     int var = ++highestVariable;
//                                     solver->add(edges[vertexOrderingColoring[i]][vertexOrderingColoring[j]]);
//                                     solver->add(-var);
//                                     solver->add(0);

//                                     solver->add(edges[vertexOrderingColoring[i]][vertexOrderingColoring[k]]);
//                                     solver->add(-var);
//                                     solver->add(0);

//                                     solver->add(edges[vertexOrderingColoring[j]][vertexOrderingColoring[k]]);
//                                     solver->add(-var);
//                                     solver->add(0);

//                                     solver->add(color1[j]);
//                                     solver->add(-var);
//                                     solver->add(0);

//                                     solver->add(color1[k]);
//                                     solver->add(-var);
//                                     solver->add(0);

//                                     solver->add(-edges[vertexOrderingColoring[i]][vertexOrderingColoring[j]]);
//                                     solver->add(-edges[vertexOrderingColoring[i]][vertexOrderingColoring[k]]);
//                                     solver->add(-edges[vertexOrderingColoring[j]][vertexOrderingColoring[k]]);
//                                     solver->add(-color1[j]);
//                                     solver->add(-color1[k]);
//                                     solver->add(var);
//                                     solver->add(0);

//                                     adjacentTriangleAndColor1.push_back(var);
//                                 }
//                             }

//                             // \lnot available1 iff or( adjacentTriangleAndColor1)
//                             for (auto l : adjacentTriangleAndColor1)
//                             {
//                                 solver->add(-l);
//                                 solver->add(-available1);
//                                 solver->add(0);
//                             }
//                             for (auto l : adjacentTriangleAndColor1)
//                                 solver->add(l);
//                             solver->add(available1);
//                             solver->add(0);

//                             // if color 0 is not available but color 1 then give color one.
//                             solver->add(available0);
//                             solver->add(-available1);
//                             solver->add(color1[i]);
//                             solver->add(0);

//                             // if no available then not color 1
//                             solver->add(+available1);
//                             solver->add(-color1[i]);
//                             solver->add(0);

//                             // if color 0 then not color 1
//                             solver->add(-color0[i]);
//                             solver->add(-color1[i]);
//                             solver->add(0);

//                             // if color then not uncolored
//                             solver->add(-color0[i]);
//                             solver->add(-uncolored[i]);
//                             solver->add(0);

//                             solver->add(-color1[i]);
//                             solver->add(-uncolored[i]);
//                             solver->add(0);
//                         }

//                         // at least one uncolored vertex
//                         for (int i = 0; i < vertices; i++)
//                             solver->add(uncolored[i]);
//                         solver->add(0);
//                     }
//                     continue;
//                 }

//                 // check 101 colorability
//                 if (non010colorable)
//                 {
//                     printf("Check colorability\n");
//                     coloring_t coloring(vertices);
//                     if (get010Coloring(vertices, matrix, coloring))
//                     {
//                         vector<lit_t> clause = get010ColoringClause(coloring, vertices, edges, triangles);
//                         for (auto lit : clause)
//                             solver->add(lit);
//                         solver->add(0);

//                         if (addPermutedColorings)
//                         {
//                             for (int i = 0; i < vertices; i++)
//                                 for (int j = i + 1; j < vertices; j++)
//                                 {
//                                     if (coloring[i] != coloring[j])
//                                     {
//                                         swap(coloring[i], coloring[j]);
//                                         vector<lit_t> clause = get010ColoringClause(coloring, vertices, edges, triangles);
//                                         for (auto lit : clause)
//                                             solver->add(lit);
//                                         solver->add(0);
//                                         swap(coloring[i], coloring[j]); // undo
//                                     }
//                                 }
//                         }
//                         stats.timeCheckFullGraphs += clock() - start;
//                         continue;
//                     }
//                 }

//                 // exclude graph from search space
//                 for (int i = 0; i < vertices; i++)
//                     for (int j = i + 1; j < vertices; j++)
//                     {
//                         if (matrix[i][j] == truth_value_true)
//                         {
//                             solver->add(-edges[i][j]);
//                         }
//                         else
//                         {
//                             solver->add(+edges[i][j]);
//                         }
//                     }
//                 solver->add(0);

//                 stats.timeCheckFullGraphs += clock() - start;

//                 printf("Solution: %d\n", ++solutionCounter);
//                 printAdjacencyMatrix(matrix);
//             }
//             else
//             {
//                 EXIT_UNWANTED_STATE
//             }
//         }
//     }
