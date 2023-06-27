#ifndef SOLVE_GENERAL_H
#define SOLVE_GENERAL_H

#include "useful.h"
#include "graphChecker.hpp"
#include "minimalityCheck.hpp"

#define COUNT_ASSIGNED 0

typedef struct
{
  int vertices;
  bool printStats;
  bool printIntermediateStatistic;
  bool hideGraphs;
  bool allModels;
  bool printFullMatrix;

  bool bipartite;
  int b_vertices[2]; // if dealing with bipartite graphs (like incidence graphs of hypergraphs), these are the partition sizes
  bool hypermode;    // for hypergraphs and hyperspace travel

  vector<int> observedVars; // all observed variables

  bool turnOffInprocessing;
  bool printFullyDefinedGraphs;

  int assignmentScoring;           // 0 = count the nubmer of assigned variables ; 1 = estimate the search space cut via recorded frequences
  int assignmentCutoff;            // create cube when a certain number of vertex pairs is decided whether the are and edge
  int assignmentCutoffPrerun;      // Number of calls of propagator before assignmentCutoff is activated
  long assignmentCutoffPrerunTime; // Seconds before calling assignmentCutoff

  bool chechSolutionInProp = true; // check solution in the propagator and continous search without restart; not applicable for more complex propagators use the incremental interface
  bool propagateLiteralsCadical;   // whether only literals should be propagated and clauses added at conflict analysis
  bool irredundantSymClauses;      // don't let the symmetry breaking clauses be part of the clause deletion policy

  bool non010colorable; // TODO get ride of as config parameter

  string proofFile;

  string cubeFile;           // name of file containing cubes
  pair<int, int> rangeCubes; // solve all cubes between rangeCubes.fist and rangeCubes.second;

  int timeout; // timeout in seconds (for each cube)

  bool hyperedgeColoring; // also excludes intersection graph

  vector<pair<int, int>> intervallsColoring; // intervalls for partitions based on coloring

  int printPartiallyDefined;

  vector<vector<lit_t>> edges;                    // edges to variables
  vector<vector<lit_t>> edges_intersection_graph; // edges of intersection graph to variables
  vector<vector<lit_t>> staticPartition;          // variables indicating whether two vertices are in the same partition
  vector<vector<vector<lit_t>>> triangles;        // variables indicating whether triangle is present

  int nextFreeVariable;

  FILE *symBreakClausesFile;
  FILE *addedClauses;
  FILE *addedColoringClauses;

  bool quiet = true; // don't print any output

  // symmetry breaking related part
#define DEFAULT_FREQUENCY 20
  int frequency = DEFAULT_FREQUENCY;
  int cutoff;                               // cutoff for minimality check
  vector<bool> initialPartition;            // initial partition for minimality check
  int frequencyConnectedComponentsSwap = 0; // based on a special coloring encoding
  bool turnoffSMS;
  bool combineStaticPlusDynamic; // start dynamic symmetry breaking with ordered partition implicitely given by some boolean variables.
  vector<vector<vertex_t>> initialVertexOrderings;
} configSolver;

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
  void solve(); // TODO suitable return value;
  GraphSolver(configSolver config)
  {
    stats.start = clock();
    this->config = config;
    this->vertices = config.vertices;
    this->edges = config.edges;
    this->nextFreeVariable = config.nextFreeVariable;

    // add dynamic symmetry breaking
    if (!config.turnoffSMS)
    {
      if (config.combineStaticPlusDynamic)
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
        auto checker = new MinimalityChecker(config.frequency, config.initialPartition, config.initialVertexOrderings, config.cutoff);
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
      config.chechSolutionInProp = false;
    }
    complexFullyDefinedGraphCheckers.push_back(checker);
  }

  vector<vector<lit_t>> edges; // edges to variables

  vector<vector<lit_t>> edge_stats;
  vector<vector<vector<lit_t>>> triangle_stats;
  statistics stats = {}; // default value initialization

  int nextFreeVariable;

private:
  vector<PartiallyDefinedGraphChecker *> partiallyDefinedGraphCheckers;
  vector<ComplexPartiallyDefinedGraphChecker *> complexPartiallyDefinedGraphCheckers;
  vector<FullyDefinedGraphChecker *> fullyDefinedGraphCheckers;
  vector<ComplexFullyDefinedGraphChecker *> complexFullyDefinedGraphCheckers;

protected:
  configSolver config;
  int vertices;

  int nModels = 0;
  const vector<int> *model = NULL; // model of the last canidate solution

  // functions which must be implemented for the concrete solver
  virtual void solve(vector<int> assumptions) = 0;              // solve the formula under the assumption
  virtual bool solve(vector<int> assumptions, int timeout) = 0; // solve with a given timeout; return false if timeout was reached

  virtual void addClause(const vector<lit_t> &clause, bool redundant) = 0;
  virtual adjacency_matrix_t getAdjacencyMatrix() = 0;
  virtual vector<truth_value_t> &getCurrentAssignemnt() = 0;
  virtual vector<vector<truth_value_t>> getStaticPartition() = 0;

  // functions which are the same for all solvers, which use the previous funcitons
private:
  bool cutoffFunction(); // If certain number of edge variables is assigned, a cube will be generated
  bool checkPartiallyDefined(bool isFullyDefined);
  bool checkFullyDefinedGraph(const adjacency_matrix_t &matrix, const vector<int> &model); // check the property of the fully defined graph, given the model
  void printStatistics();

  void recordGraphStats(const adjacency_matrix_t &matrix);
  void initEdgeMemory();
  void initTriangleMemory();
  void printEdgeStats();

public:
  bool propagate(); // Check state of partial assignment and add clauses if necessary; returns true if no clause was added otherwise false
  bool check();     // prints and excludes graph from search space; returns true if no clause was added otherwise false
};

#endif
