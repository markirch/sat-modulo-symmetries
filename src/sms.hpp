#ifndef SOLVE_GENERAL_H
#define SOLVE_GENERAL_H

#include <cassert>

#include "useful.h"
#include "graphChecker.hpp"
#include "minimalityCheck.hpp"

using std::pair;
using std::string;

enum
{
  COUNT_ASSIGNED = 0,
  COUNT_ASSIGNED_WEIGHTED = 1,
};

class SolverConfig
{
public:
  int vertices = 2;
  bool printStats = false;
  bool printIntermediateStatistic = false;
  bool hideGraphs = false;
  bool allModels = false;
  bool printFullMatrix = false;
  bool printFullModel = false;

  bool bipartite = false;
  int b_vertices[2];      // if dealing with bipartite graphs (like incidence graphs of hypergraphs), these are the partition sizes
  bool hypermode = false; // for hypergraphs and hyperspace travel

  vector<int> observedVars; // all observed variables

  bool turnOffInprocessing = false;
  bool printFullyDefinedGraphs = false;
  int printPartiallyDefined = 0;

  int assignmentScoring = COUNT_ASSIGNED;
  int assignmentCutoff = 0;            // create a cube when the assignment score exceeds this
  int assignmentCutoffPrerun = 0;      // number of calls of propagator before assignmentCutoff is activated
  long assignmentCutoffPrerunTime = 0; // time in seconds before assignmentCutoff is activated

  bool checkSolutionInProp = true;       // check solution in the propagator and continue search without restarting; for more complex propagators, use the incremental interface
  bool propagateLiteralsCadical = false; // whether only literals should be propagated and clauses added at conflict analysis
  bool irredundantSymClauses = false;    // prevent symmetry breaking clauses from being deleted on clause DB cleanup

  bool non010colorable = false; // TODO get rid of as config parameter

  string proofFile;

  string cubeFile;           // name of file containing cubes
  pair<int, int> rangeCubes; // solve all cubes between rangeCubes.fist and rangeCubes.second;

  int timeout = 0; // timeout in seconds (for each cube)

  bool hyperedgeColoring = false; // also excludes intersection graph

  vector<pair<int, int>> intervallsColoring; // intervalls for partitions based on coloring

  vector<vector<lit_t>> edges;                    // edges to variables
  vector<vector<lit_t>> edges_intersection_graph; // edges of intersection graph to variables
  vector<vector<lit_t>> staticPartition;          // variables indicating whether two vertices are in the same partition
  vector<vector<vector<lit_t>>> triangles;        // variables indicating whether triangle is present

  int numberOfOverlayingGraphs = 0;            // generate multiple graphs at once (can also bee seen as edge coloring given by the different level)
  vector<vector<vector<lit_t>>> edgesMultiple; // edge variables of several graphs

  int nextFreeVariable = 1;

  FILE *symBreakClausesFile = NULL;
  FILE *addedClauses = NULL;
  FILE *addedColoringClauses = NULL;

  bool quiet = true; // don't print any output

  // symmetry breaking related part
#define DEFAULT_FREQUENCY 20
  int frequency = DEFAULT_FREQUENCY;
  int cutoff = 0;                           // cutoff for minimality check
  vector<bool> initialPartition;            // initial partition for minimality check
  int frequencyConnectedComponentsSwap = 0; // based on a special coloring encoding
  bool turnoffSMS = false;
  bool combineStaticPlusDynamic = false; // start dynamic symmetry breaking with ordered partition implicitely given by some boolean variables.
  vector<vector<vertex_t>> initialVertexOrderings;

  // DEFAULT CONSTRUCTORS

  // without any arguments, construct a dummy config that turns off SMS
  // and essentially acts like an ordinary SAT solver
  // WARNING: there must be at least 2 vertices, otherwise there's no edge and a segfault

  SolverConfig(int vertices = 2)
  {
    assert(vertices >= 2);
    set_vertices(vertices);
    if (vertices == 2)
    {
      turnoffSMS = true;
    }
  }

  SolverConfig(int vertices, int cutoff) : cutoff(cutoff) {
    set_vertices(vertices);
  }

  // change the number of vertices in an existing SolverConfig object this way to update all dependencies
  void set_vertices(int vertices);

  void init_edge_vars();
  void init_multi_edge_vars(); // to be called only after init_edge_vars()
  void init_triangle_vars(int triangleVars = 0);
  void init_intersection_vars(int &minIntersectionVar);
};

void make_multi_edge_vars(SolverConfig &config);

typedef struct
{
  clock_t start;
  clock_t timePropagator;
  clock_t timeCheckFullGraphs;
  long long callsPropagator;
  long long callsCheck;
  long long hyperclauses;
} statistics;

class GraphSolver
{
public:
  bool solve(); // TODO suitable return value;
  GraphSolver(SolverConfig &config)
  {
    stats.start = clock();
    this->config = config;
    this->vertices = config.vertices;
    this->edges = config.edges;
    this->nextFreeVariable = config.nextFreeVariable;

    // add dynamic symmetry breaking
    if (!config.turnoffSMS)
    {
      if (config.numberOfOverlayingGraphs)
      {
        if (config.initialVertexOrderings.empty())
        {
          vertex_ordering_t basicVertexOrdering;
          for (int i = 0; i < vertices; i++)
            basicVertexOrdering.push_back(i);
          config.initialVertexOrderings.push_back(basicVertexOrdering);
        }

        auto checker = new MultipleMinimalityChecker(config.frequency, config.initialPartition, config.initialVertexOrderings, config.cutoff, config.symBreakClausesFile);
        this->partiallyDefinedMultiGraphCheckers.push_back(checker);
      }
      else if (config.combineStaticPlusDynamic)
      {
        auto checker = new MinimalityCheckerWithStaticPartition(config.frequency, config.cutoff, config.edges, config.staticPartition);
        this->complexPartiallyDefinedGraphCheckers.push_back(checker);
      }
      else
      {

        if (config.initialVertexOrderings.empty())
        {
          vertex_ordering_t basicVertexOrdering;
          for (int i = 0; i < vertices; i++)
            basicVertexOrdering.push_back(i);
          config.initialVertexOrderings.push_back(basicVertexOrdering);
        }
        auto checker = new MinimalityChecker(config.frequency, config.initialPartition, config.initialVertexOrderings, config.cutoff, config.symBreakClausesFile);
        this->partiallyDefinedGraphCheckers.push_back(checker);
      }
    }
  }
  virtual ~GraphSolver() {}

  // additional properties which are ensured thru some propagators

  void addPartiallyDefinedGraphChecker(PartiallyDefinedGraphChecker *checker)
  {
    partiallyDefinedGraphCheckers.push_back(checker);
  }

  void addComplexPartiallyDefinedGraphChecker(ComplexPartiallyDefinedGraphChecker *checker)
  {
    complexPartiallyDefinedGraphCheckers.push_back(checker);
  }

  void addFullyDefinedGraphChecker(FullyDefinedGraphChecker *checker)
  {
    fullyDefinedGraphCheckers.push_back(checker);
  }

  void addComplexFullyDefinedGraphChecker(ComplexFullyDefinedGraphChecker *checker)
  {
    if (!checker->addsOnlyObservedLiterals)
    {
      printf("Warning: uses incremental interface for adding clauses from special propagator\n");
      config.checkSolutionInProp = false;
    }
    complexFullyDefinedGraphCheckers.push_back(checker);
  }

  void addPartiallyDefinedMultiGraphChecker(PartiallyDefinedMultiGraphChecker *checker)
  {
    partiallyDefinedMultiGraphCheckers.push_back(checker);
  }

  clause_t theClauseThatBlocks(const forbidden_graph_t &);

  vector<vector<lit_t>> edges; // edges to variables

  vector<vector<lit_t>> edge_stats;
  vector<vector<vector<lit_t>>> triangle_stats;
  statistics stats = {}; // default value initialization

  int nextFreeVariable;
  SolverConfig config;

private:
  vector<PartiallyDefinedGraphChecker *> partiallyDefinedGraphCheckers;
  vector<ComplexPartiallyDefinedGraphChecker *> complexPartiallyDefinedGraphCheckers;
  vector<FullyDefinedGraphChecker *> fullyDefinedGraphCheckers;
  vector<ComplexFullyDefinedGraphChecker *> complexFullyDefinedGraphCheckers;
  vector<PartiallyDefinedMultiGraphChecker *> partiallyDefinedMultiGraphCheckers;

protected:
  int vertices;

  int nModels = 0;
  const vector<int> *model = NULL; // model of the last canidate solution

  // functions which must be implemented for the concrete solver
  virtual bool solve(vector<int> assumptions) = 0;              // solve the formula under the assumption
  virtual bool solve(vector<int> assumptions, int timeout) = 0; // solve with a given timeout; return false if timeout was reached

  virtual adjacency_matrix_t getAdjacencyMatrix() = 0;
  virtual vector<adjacency_matrix_t> getAdjacencyMatrixMultiple() = 0;
  virtual vector<truth_value_t> &getCurrentAssignemnt() = 0;
  virtual vector<vector<truth_value_t>> getStaticPartition() = 0;

  virtual void printFullModel() = 0;

  void recordGraphStats(const adjacency_matrix_t &matrix);
  void initEdgeMemory();
  void initTriangleMemory();
  void printEdgeStats();

  // functions which are the same for all solvers, which use the previous funcitons
private:
  bool cutoffFunction(); // If certain number of edge variables is assigned, a cube will be generated
  bool checkPartiallyDefined(bool isFullyDefined);
  bool checkFullyDefinedGraph(const adjacency_matrix_t &matrix, const vector<int> &model); // check the property of the fully defined graph, given the model

public:
  bool propagate(); // Check state of partial assignment and add clauses if necessary; returns true if no clause was added otherwise false
  bool check();     // prints and excludes graph from search space; returns true if no clause was added otherwise false
  virtual void addClause(const vector<lit_t> &clause, bool redundant) = 0;
  void printStatistics();
};

#endif
