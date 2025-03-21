/**
 * Implementation of common parts amongst all solvers
 */
#include <cmath>
#include "sms.hpp"
#include "minimalityCheck.hpp"

GraphSolver::GraphSolver(SolverConfig config, struct minimality_config_t minimality_config) : config(config)
{
  current_trail.push_back(vector<int>());
  this->vertices = config.vertices;
  connect_external_propagator(this);
  connect_fixed_listener(this);
  graphHandler = new GraphHandler(vertices, config.directed);
  int numVariables = graphHandler->getNumVariables();
  this->numVars = numVariables;
  currentAssignment = vector<truth_value_t>(numVariables + 1, truth_value_unknown);
  isFixed = vector<bool>(numVariables + 1, false);

  // this->set("binary", 0);
  // this->set("lidrup", 1);
  // this->trace_proof("qwer");
  // this->set("flushproof", 1);

  if (config.assignmentCutoff || config.createGame)
  {
    this->set("ilb", 1);
    this->set("ilbassumptions", 1);
  }

  for (int i = 1; i <= numVariables; i++)
    add_observed_var(i);

  if (!minimality_config.turnoffSMS)
  {
    // add minimality checker
    addPartiallyDefinedGraphChecker(new MinimalityChecker(minimality_config, vertices, config.directed));
  }

  // TODO adapt on whether directed/bipartite/hypergraph ....

  // if (config.createGame || config.assignmentCutoff)
  // {
  //   if (!this->set("ilb", 1))
  //     EXIT_UNWANTED_STATE
  //   if (!this->set("ilbassumptions", 1))
  //     EXIT_UNWANTED_STATE
  // }

  if (!config.cadicalConfig.empty())
  {
    for (auto option : config.cadicalConfig)
    {
      std::cout << "Setting Cadical options '" << option << "'" << std::endl;
      // if option is not starting with `--` then add it
      if (option.substr(0, 2) != "--")
        option = "--" + option;
      if (!this->set_long_option(option.c_str()))
      {
        std::cerr << "invalid Cadical option '" << option << "'" << std::endl;
        throw std::runtime_error("Invalid Cadical option");
      }
    }
  }
}

// TODO be careful with fully defined graph checkers to adapt observed and increase arrays if necessary

bool GraphSolver::checkPartiallyDefined(bool isFullyDefined)
{
  if (inLookaheadState)
    return true;

  adjacency_matrix_t matrix = graphHandler->assignment2graph(currentAssignment);
  try
  {
    PRINT_CURRENT_LINE
    for (auto checker : this->partiallyDefinedGraphCheckers)
      checker->check(matrix, isFullyDefined);
  }
  catch (const forbidden_graph_t forbiddenGraph)
  {
    addClause(graphHandler->theClauseThatBlocks(forbiddenGraph), config.redundantPropClauses);
    return false;
  }

  if (!inPrerunState && config.simpleAssignmentCutoff && graphHandler->numAssigned(currentAssignment) >= config.simpleAssignmentCutoff) // TODO also check that not in prerun state
  {
    addClause(graphHandler->solutionBlockingClause(currentAssignment), false);
    graphHandler->printCube(currentAssignment);
    return false;
  }

  return true;
}

// TODO maybe only declare all the class here and give the implementation at the end, because not that important.

class TimeoutTerminator : public CaDiCaL::Terminator
{
private:
  int timeout;
  clock_t start_time;

public:
  TimeoutTerminator(int timeout) : timeout(timeout)
  {
    start_time = clock();
  }

  bool terminate() override
  {
    return ((double)clock() - start_time) / CLOCKS_PER_SEC > timeout;
  }
};

// Allows parsing a cube file and returns the next cube (in the given range)
class CubeReader
{
private:
  std::ifstream cubeFile;
  std::pair<int, int> cubeRange;
  int currentLine = 1;

public:
  CubeReader(const std::string &filename, const vector<int> range, int cubeLine)
  {
    cubeFile.open(filename);
    if (!cubeFile.is_open())
    {
      throw std::runtime_error("Could not open cube file");
    }

    if (cubeLine != 0)
    {
      if (!range.empty())
        EXIT_UNWANTED_STATE
      cubeRange = {cubeLine, cubeLine};
    }
    else if (range.size() == 2)
      cubeRange = {range[0], range[1]};
    else
    {
      throw std::runtime_error("Invalid range for cubes");
    }

    // ignore all lines until the first cube
    std::string line;
    while (currentLine != cubeRange.first)
    {
      currentLine++;
      std::getline(cubeFile, line);
    }
  }

  ~CubeReader()
  {
    if (cubeFile.is_open())
    {
      cubeFile.close();
    }
  }

  bool hasNextCube()
  {
    return currentLine <= cubeRange.second && !cubeFile.eof();
  }

  std::vector<int> nextCube(int &cubeNr)
  {
    cubeNr = currentLine;
    std::string line;
    std::vector<int> cube;
    std::getline(cubeFile, line);
    currentLine++;
    std::istringstream iss(line);
    std::string lit;
    while (iss >> lit)
    {
      if (lit == "a" || lit == "0")
      {
        continue;
      }
      cube.push_back(std::stoi(lit));
    }

    return cube;
  }
};

int GraphSolver::sms_solve()
{
  PRINT_CURRENT_LINE
  stats.start = clock();

  if (config.prerunTime)
  {
    inPrerunState = true;
    vector<int> assumptions;
    int res = sms_main_loop(assumptions, config.prerunTime);
    inPrerunState = false;

    if (res != 0) // return value 0 means unknown
    {
      LOG(LOG_LEVEL_INFO, "Instance already solved during prerun");
      return res;
    }
  }

  if (!config.simplifiedFormulaFile.empty() || !config.learnedClausesFile.empty())
    this->simplify(); // do some preprocessing

  if (!config.simplifiedFormulaFile.empty())
  {
    writeSimplified(config.simplifiedFormulaFile);

    if (config.learnedClausesFile.empty())
      return 0;
  }

  if (!config.learnedClausesFile.empty())
  {
    writeLearnedClauses(config.learnedClausesFile);
    return 0;
  }

  if (config.assignmentCutoff)
  {
    generateCubes(config.assignmentCutoff);
    return 0;
  }

  if (config.createGame)
  {
    vector<int> assumptions;
    if (!config.cubeFile.empty())
    {
      CubeReader r(config.cubeFile, config.cubesRange, config.cubeLine);

      while (r.hasNextCube())
      {
        int cubeNr;
        vector<int> assumptions = r.nextCube(cubeNr);
        createGame(config.createGameProp, assumptions, config.createGameRecLvl);
      }
    }
    else
    {
      createGame(config.createGameProp, vector<int>(), config.createGameRecLvl);
    }
    return 0;
  }

  if (!config.cubeFile.empty()) // check whether it should be solved with cubes as assumptions
  {
    CubeReader r(config.cubeFile, config.cubesRange, config.cubeLine);

    while (r.hasNextCube())
    {

      int cubeNr;
      vector<int> assumptions = r.nextCube(cubeNr);
      printf("Solve cube %d\n", cubeNr);
      // printf("Cube:");
      // for (int lit : assumptions)
      //   printf(" %d", lit);
      // printf("\n");

      clock_t start = clock();
      if (sms_main_loop(assumptions, config.cubeTimeout) == 0)
      {
        printf("Timeout reached for solving cube\n");
      }
      else
      {
        printf("Time for cube %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
      }
      printStats();
    }

    printf("All cubes processed\n");
    return 0;
  }

  vector<int> assumptions;
  int res = sms_main_loop(assumptions, config.timeout);
  if (res == 0)
  {
    LOG(LOG_LEVEL_INFO, "Instance is unknown");
  }
  printStats();

  LOG(LOG_LEVEL_INFO, "Result: " << res);
  return res;
}

int GraphSolver::sms_main_loop(vector<int> assumptions, int timeout)
{
  PRINT_CURRENT_LINE
  TimeoutTerminator *terminator;
  if (timeout)
  {
    terminator = new TimeoutTerminator(timeout);
    this->connect_terminator(terminator);
  }

  int res;
  do
  {
    // assumptions have to be added before each call of solve
    for (int lit : assumptions)
      CaDiCaL::Solver::assume(lit);
    res = CaDiCaL::Solver::solve();
    if (res == 10)
    {
      // check complex fully defined graph checkers
      if (!checkIncremental())
        continue;

      if (!config.hideGraphs)
      {
        adjacency_matrix_t matrix = graphHandler->assignment2graph(currentAssignment); // TODO check whether crrentAssignment can really be used outside IPASIR-UP
        graphHandler->print(matrix);
      }

      if (config.allModels)
        blockSolution();
      else
        break;
    }
  } while (res == 10);

  if (timeout)
  {
    this->disconnect_terminator();
    delete terminator;
  }
  return res;
}

bool GraphSolver::checkIncremental()
{
  vector<int> model = {0};
  for (int i = 1; i <= this->vars(); i++)
    model.push_back(this->val(i) > 0 ? i : -i);

  if (this->vars() != this->numVars)
  {
    LOG(LOG_LEVEL_INFO, "Number of variables in cadical doesn't coinside with the number of variables of the graph solver. All remaining variables are assumed to be true");
    LOG(LOG_LEVEL_INFO, "Number of variables in cadical: " + std::to_string(this->vars()));
    LOG(LOG_LEVEL_INFO, "Number of variables in graph solver: " + std::to_string(this->numVars));
    for (int i = this->vars() + 1; i <= this->numVars; i++)
    {
      model.push_back(i);
    }
  }

  bool res = true;

  int nextFreeVariable = this->numVars + 1;
  adjacency_matrix_t matrix = graphHandler->assignment2graph(currentAssignment);
  try
  {
    for (auto checker : this->complexFullyDefinedGraphCheckers)
      checker->check(matrix, model, nextFreeVariable);
  }
  catch (const vector<clause_t> clauses)
  {
    for (clause_t c : clauses)
      addClause(c, false);

    res = false;
  }

  this->numVars = std::max(this->numVars, nextFreeVariable - 1);

  return res;
}

void GraphSolver::blockSolution()
{
  adjacency_matrix_t matrix = graphHandler->assignment2graph(currentAssignment);
  stats.nModels++;

  vector<lit_t> clause = graphHandler->solutionBlockingClause(currentAssignment);
  addClause(clause, false);
}

void GraphSolver::addClause(const vector<lit_t> &clause, bool redundant)
{
  PRINT_CURRENT_LINE

  CaDiCaL::State s = state();
  int isReady = s & CaDiCaL::State::READY; // means that it is not currently solving and clauses can be added normally. TODO check again
  if (!isReady)
  {
    clauses.push_back(make_pair(clause, redundant));
  }
  else
  {
    PRINT_CURRENT_LINE
    // use incremental interface
    for (auto l : clause)
      this->add(l);
    this->add(0);
  }
}

//  ------------ implementation of IPASIR-UP interface ------------

void GraphSolver::notify_assignment(const std::vector<int> &lits)
{
  for (auto lit : lits)
  {
    changeInTrail = true;
    int absLit = abs(lit);
    currentAssignment[absLit] = lit > 0 ? truth_value_true : truth_value_false;
    // this->isFixed[absLit] = is_fixed;
    current_trail.back().push_back(lit);
  }
}

void GraphSolver::notify_new_decision_level()
{
  current_trail.push_back(vector<int>());

  stats.maxDepth = std::max(stats.maxDepth, (int)current_trail.size());
  stats.summedDepth += current_trail.size();
  stats.depthCalls++;
}

void GraphSolver::notify_backtrack(size_t new_level)
{
  while (current_trail.size() > new_level + 1)
  {
    auto last = current_trail.back();
    for (int l : last)
    {
      if (!isFixed[abs(l)])
        currentAssignment[abs(l)] = truth_value_unknown;
    }
    current_trail.pop_back();
  }
}

// currently not checked in propagator but with the normal incremental interface to allow adding other literals or even new once.
bool GraphSolver::cb_check_found_model(const std::vector<int> &)
{
  PRINT_CURRENT_LINE
  if (!clauses.empty())
    return false;

  if (!checkPartiallyDefined(true))
    return false;

  if (config.allModels && complexFullyDefinedGraphCheckers.empty())
  {
    PRINT_CURRENT_LINE
    if (!config.hideGraphs)
    {
      adjacency_matrix_t matrix = graphHandler->assignment2graph(currentAssignment);
      graphHandler->print(matrix);
    }

    blockSolution();
    return false;
  }

  return true;
}

bool GraphSolver::cb_has_external_clause(bool &redundant)
{
  // printf("ASDF\n");
  // printf("Trail %ld:", current_trail.size());
  // for (auto lits : current_trail)
  // {
  //   for (auto lit : lits)
  //     printf(" %d", lit);
  //   printf(" ; ");
  // }
  // printf("\n");

  if (clauses.empty() && changeInTrail)
  {
    changeInTrail = false;
    checkPartiallyDefined(false);
  }

  if (!clauses.empty())
  {
    redundant = clauses.back().second;
    return true;
  }
  return false;
}

int GraphSolver::cb_add_external_clause_lit()
{
  vector<int> &lastClause = clauses.back().first;
  if (lastClause.empty())
  {
    clauses.pop_back(); // delete last clause
    return 0;
  }
  else
  {
    int lit = lastClause.back();
    lastClause.pop_back();
    return lit;
  }
}

// --------------- statistics ---------------

void GraphSolver::printStats()
{
  if (config.allModels)
    printf("Number of graphs: %d\n", stats.nModels);

  printf("Maximal depth: %d\n", stats.maxDepth);
  printf("Average depth: %f\n", (double)stats.summedDepth / stats.depthCalls);

  for (auto checker : this->partiallyDefinedGraphCheckers)
    checker->printStats();

  for (auto checker : this->complexFullyDefinedGraphCheckers)
    checker->printStats();
}

// ------------------other------------------------

// TODO maybe merge ExtendedClauseCounter and ExtendedClauseWriter

class ExtendedClauseCounter : public CaDiCaL::ClauseIterator
{
public:
  int count = 0;
  int countRedundant = 0;
  int countRedundantSmall = 0;
  int maxVar = 0;

  int maxRedundantClausesize = 0;

  ExtendedClauseCounter(int maxRedundantClausesize = 0) : maxRedundantClausesize(maxRedundantClausesize) {}

  bool redundant_clause(const std::vector<int> &c) override
  {
    countRedundant++;
    if ((int)c.size() > maxRedundantClausesize)
      return true; // skip too large learned clauses
    countRedundantSmall++;
    return clause(c);
  }

  bool clause(const std::vector<int> &c) override
  {
    count++;
    for (auto l : c)
      maxVar = std::max(maxVar, abs(l));
    return true;
  }
};

class ExtendedClauseWriter : public CaDiCaL::ClauseIterator
{
public:
  FILE *file;
  int maxRedundantClausesize = 0;

  ExtendedClauseWriter(FILE *file, int maxRedundantClausesize = 0) : file(file), maxRedundantClausesize(maxRedundantClausesize) {};

  bool redundant_clause(const std::vector<int> &c) override
  {
    if ((int)c.size() > maxRedundantClausesize)
      return true; // skip too large learned clauses
    return clause(c);
  };

  bool clause(const std::vector<int> &c) override
  {
    for (auto l : c)
      fprintf(file, "%d ", l);
    fprintf(file, "0\n");
    return true;
  };
};

class LearnedClauseWriter : public CaDiCaL::ClauseIterator
{
public:
  FILE *file;
  int maxRedundantClausesize = 5; // TODO make adaptable later

  LearnedClauseWriter(FILE *file) : file(file) {};

  bool redundant_clause(const std::vector<int> &c) override
  {
    if ((int)c.size() > maxRedundantClausesize)
      return true; // skip too large learned clauses

    for (auto l : c)
      fprintf(file, "%d ", l);
    fprintf(file, "0\n");
    return true;
  };

  bool clause(const std::vector<int> &) override
  {
    return true; // nothing to do
  };
};

void GraphSolver::writeSimplified(const string filename)
{
  ExtendedClauseCounter counter(config.maxPrintedLearnedClauseSize);
  traverse_clauses(counter, true, true);

  FILE *file = fopen(filename.c_str(), "w");
  if (!file)
    EXIT_UNWANTED_STATE
  // print p cnf number_of_variables number_of_clauses
  fprintf(file, "p cnf %d %d\n", counter.maxVar, counter.count);
  ExtendedClauseWriter writer(file, config.maxPrintedLearnedClauseSize);
  traverse_clauses(writer, true, true);
  fclose(file);
}

void GraphSolver::writeLearnedClauses(const string filename)
{
  FILE *file = fopen(filename.c_str(), "w");
  if (!file)
    EXIT_UNWANTED_STATE
  // print p cnf number_of_variables number_of_clauses
  LearnedClauseWriter writer(file);
  traverse_clauses(writer, true, true);
  fclose(file);
}