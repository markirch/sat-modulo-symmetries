#ifndef SOLVE_GENERAL_H
#define SOLVE_GENERAL_H

#include <cassert>

#include "useful.h"
#include "graph.hpp"

#include "graphChecker.hpp"
#include "minimalityCheck.hpp"
#include "cadical.hpp"

using std::pair;
using std::string;

typedef struct
{
  // graph structure
  int vertices = 1;
  bool directed = false;
  bool allModels = false;
  bool hideGraphs = false; // if true, the graphs are not printed
  bool checkSolutionInProp = true;
  bool redundantPropClauses = false; // if true clauses added by propagators are treated as redundant
  int timeout = 0;                   // timeout in seconds

  int prerunTime = 0;           // time in seconds to run before creating cubes/output simplify formula/ output learned clauses/ activating lookahead.
  string simplifiedFormulaFile; // write simplified formula to the file including small redundant clauses TODO be more specific
  string learnedClausesFile;    // write all learned clauses up to the given size to the file TODO check whether this is also for irredundant clauses
  int maxPrintedLearnedClauseSize = 5;

  // cubing related
  int assignmentCutoff = 0;       // create a cube when the number of assigned variables exceeds this
  int simpleAssignmentCutoff = 0; // create cubes soley based on edge variables. Cubes are blocked by clauses until unsat is reached
  string cubeFile;                // filename containing the cubes
  string cubeFileTest;            // test the given cubes, i.e., whether negation of the cubes is unsat
  int cubeLine = 0;               // the cube from the given line is selected (starting with line 1)
  vector<int> cubesRange;         // range of cubes to be processed
  int cubeTimeout = 0;
  bool cubesOnlyDecisions = false; // the output of the cubes only contains the decision literals and not implied literals
  bool lookaheadOnlyEdgeVars = false;

  bool createGame = false;  // print the trail on a conflict or if all clauses are satisfied using lookeahead. Afterwards, the solving stops
  float createGameProp = 0; // randomize the truth value of the branching variable. The value is the probability of flipping the truth value.
  int createGameRecLvl = 0; // the depth of the lookahead on the truth values for creating games

  vector<string> cadicalConfig; // space-separated list of command-line options for cadical without the -- prefixes
} SolverConfig;

typedef struct
{
  clock_t start;
  clock_t timePropagator;
  clock_t timeCheckFullGraphs;
  long long callsPropagator;
  long long callsCheck;
  int nModels;
  int maxDepth = 0;
  long long summedDepth = 0; // sum of the depth of all calls of propagate
  long long depthCalls = 0;
} sms_statistics_t;

class GraphSolver : public CaDiCaL::ExternalPropagator, public CaDiCaL::FixedAssignmentListener, public CaDiCaL::Solver
{
public:
  GraphSolver(SolverConfig config, struct minimality_config_t minimality_config);
  ~GraphSolver()
  {
    if (graphHandler)
      delete graphHandler;
    for (auto c : partiallyDefinedGraphCheckers)
      delete c;
    for (auto c : complexFullyDefinedGraphCheckers)
      delete c;
  }

  SolverConfig config;
  int numVars = 0;

  int sms_solve();                                                                     // TODO Same return values as CaDiCaL::solve but also takes fully defined graphs into account
  vector<vector<int>> getEdgeVariables() { return graphHandler->getEdgeVariables(); }; // returns the edge variables used for representing the graph
  void printStats();
  void generateCubes(int assignmentCutoff);
  void createGame(float randomize = 0, vector<int> assumptions = {}, int recLvl = 0); // print the trail on a conflict or if all clauses are satisfied using lookeahead

private:
  int countAssignedOrInactive(const vector<vector<int>> &cubeTrail); // count the number of assigned variables and inactive variables
  
  typedef struct
  {
    vector<bool> isAssigned;       // indicate whether a variable is assigned
    vector<vector<int>> cubeTrail; // trail of assigned literals
    vector<vector<int>> decisions; // trail of variables set leading to the result (allows setting multiple decisions per level and assumptions and root level)

    bool noLearn = false; // if true then it is assumed that the solver does not learn clauses by conflict analyis
  } cubing_state_t;

  void generateCubesRec(cubing_state_t &cubingState, int assignmentCutoff);
  bool checkCutoff(const vector<vector<int>> &cubeTrail, int assignmentCutoff); // check if the number of assigned and inactive variables exceeds the cutoff
  /**
   *  For each eligable variable compute the score; if any literal fails then restart the lookahead and add negation of the literal to the cube trail.
   *
   * @param numPropagted stores the number of propagated literals for each truth assignemnt of each variable
   *
   * @return false if the current assignment cannot lead to a solution
   */
  bool lookahead(cubing_state_t &cubingState, const bool repeatOnConflict, vector<std::pair<int, int>> &numPropagted);
  bool createGameRec(cubing_state_t &cubingState, int recDepth, int &branchingLit, int &finalScore); // returns false if a final state was reached, otherwise branching literal and score are stored

  /**
   * @brief Updated the cube trail with the given literals including propagated literals. If 20 is returned then the cube trail is not updated.
   *
   * @param lits A list of lits which are added to the cube trail
   * @param isAssigned Indicates whether a variable is already assigned
   * @param cubeTrail The current state, i.e., trail, assigned variables, and decisions
   * @param numPropagated If the propagation was successful then the number of propagations is stored here (excluding fixed assignments). If UNSAT then the number of all variables
   * @return int The result of the propagation (0 if UNKNOWN, 10 if SAT, 20 if UNSAT)
   */
  int updateCubeTrail(const vector<int> lits, cubing_state_t &cubingState, int &numPropagated);
  void undoCubeTrail(cubing_state_t &cubingState);

private:
  int solve() { EXIT_UNWANTED_STATE } // deactivates the original solve function
  void assume(int) { EXIT_UNWANTED_STATE }

  // part for adding propagators
public:
  void addPartiallyDefinedGraphChecker(PartiallyDefinedGraphChecker *checker) { partiallyDefinedGraphCheckers.push_back(checker); }
  void addComplexFullyDefinedGraphChecker(ComplexFullyDefinedGraphChecker *checker) { complexFullyDefinedGraphCheckers.push_back(checker); }

private:
  GraphHandler *graphHandler;
  vector<PartiallyDefinedGraphChecker *> partiallyDefinedGraphCheckers;
  vector<ComplexFullyDefinedGraphChecker *> complexFullyDefinedGraphCheckers;

  // ------------------------------------propagators and listeners------------------------------------

private:
  int vertices = 1;
  sms_statistics_t stats = {};
  vector<truth_value_t> currentAssignment; // currentAssignment[v] gives the truthvalue of variable v (if observed)

  bool inLookaheadState = false; // if true then skip propagator calls
  bool inPrerunState = false;

  // additional properties which are ensured by some propagators

  int sms_main_loop(vector<int> assumptions, int timeout);

  bool checkPartiallyDefined(bool isFullyDefined); // check partially defined graph and add clauses if necessary; returns true if no clause was added otherwise false
  bool checkIncremental();                         // check complex fully defined graph checkers outside of the propagator
  void addClause(const vector<lit_t> &clause, bool redundant);

  void blockSolution(); // print combinatorial object and block the current solution by adding a clause

  // ------------ extended solver interface -------------------------------

  void writeSimplified(const string fileName);     // write simplified formula to the given file (based on the current solver state)
  void writeLearnedClauses(const string fileName); // write all learned clauses up to the given size to the file

  // ------------ implementation of fixed assignment listener ------------

private:
  vector<int> fixedLiterals;
  vector<bool> isFixed;

public:
  void notify_fixed_assignment(int lit)
  {
    // printf("Fixed assignment: %d\n", lit);
    if (abs(lit) >= isFixed.size())
    {
      isFixed.resize(abs(lit) + 1, false);
    }

    if (!isFixed[abs(lit)]) // avoid adding multiple times
    {
      isFixed[abs(lit)] = true;
      fixedLiterals.push_back(lit);
    }
  }

  // ------------ implementation of IPASIR-UP interface ------------

private:
  vector<vector<int>> current_trail;
  bool changeInTrail = true; // true if since last call of propagate the trail of the observed variables has changed

  vector<pair<vector<int>, bool>> clauses; // all clauses which should be added. The second value indicates whether the clause is treated as redundant or not

public:
  void notify_assignment(const std::vector<int> &lits);
  void notify_new_decision_level();
  void notify_backtrack(size_t new_level);

  // currently not checked in propagator but with the normal incremental interface to allow adding other literals or even new once.
  bool cb_check_found_model(const std::vector<int> &model);
  bool cb_has_external_clause(bool &redundant);
  int cb_add_external_clause_lit();

  // unused functions
  int cb_decide() { return 0; }
  int cb_propagate() { return 0; }
  int cb_add_reason_clause_lit(int) { return 0; };
};

#endif
