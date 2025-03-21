/**
 *  Implementation for creating cubes using the 'propagate' function from the CaDiCaL library.
 *  Might be slower than using the IPASIR-UP interface but much cleaner and easier to adapt.
 */

#include "useful.h"
#include "sms.hpp"
#include <set>

#define IS_FIXED(lit) (abs(lit) < (int)isFixed.size() && isFixed[abs(lit)])

void printCube(const vector<vector<int>> &cubeTrail, const vector<int> &fixedLiterals, bool onlyDecisions = true)
{
    if (!onlyDecisions)
    {
        printf("a");
        for (auto lit : fixedLiterals)
            printf(" %d", lit);
        for (auto lits : cubeTrail)
            for (auto lit : lits)
                printf(" %d", lit);
        printf(" 0\n");

        // print ordered cubes
        // vector<int> allLits;
        // // Collect all literals
        // for (auto lit : fixedLiterals)
        //     allLits.push_back(lit);
        // for (auto lits : cubeTrail)
        //     for (auto lit : lits)
        //         allLits.push_back(lit);

        // // Sort literals by absolute value
        // std::sort(allLits.begin(), allLits.end(), [](int a, int b)
        //           { return abs(a) < abs(b); });

        // // Print sorted literals
        // printf("a");
        // for (auto lit : allLits)
        //     printf(" %d", lit);
        // printf(" 0\n");
    }
    else
    {
        printf("Decisions: ");
        bool first = true;
        for (auto lits : cubeTrail)
        {
            if (first)
            {
                first = false;
                continue;
            }
            printf(" %d", lits[0]);
        }
        printf("\n");
    }

    // count without duplicates
    // std::set<int> cubeSet;
    // for (auto lits : cubeTrail)
    //     for (auto lit : lits)
    //         cubeSet.insert(lit);
    // for (auto lit : fixedLiterals)
    //     cubeSet.insert(lit);
    // printf("c Cube size: %d\n", (int)cubeSet.size());
}

void GraphSolver::generateCubes(int assignmentCutoff)
{
    this->simplify();

    // this->set("reduce", 0); // do not remove any clauses

    cubing_state_t cubingState = {};
    cubingState.isAssigned = vector<bool>(this->vars() + 1, false);

    // printf("Number of vars: %d\n", this->vars());

    int numPropagated;
    int res = updateCubeTrail(vector<int>(), cubingState, numPropagated); // propagate empty assumptions to get initial assignments
    if (res == 20)                                                        // UNSAT
    {
        LOG(LOG_LEVEL_INFO, "Instance is UNSAT");
        exit(1);
    }

    printf("c Number of variables: %d\n", vars());
    printf("c Number of active variables: %d\n", active());

    // printf("c Initial trail: ");
    // for (auto lit : cubeTrail.back())
    //     printf("%d ", lit);
    // printf("\n");

    generateCubesRec(cubingState, assignmentCutoff);
}

inline int getBestVar(const vector<std::pair<int, int>> &scores, std::pair<int, int> &varScore)
{
    double bestScore = -1;
    int bestVar = 0;
    for (int i = 1; i < (int)scores.size(); i++)
    {
        if (scores[i].first == 0)
            continue;
        double score = std::min(scores[i].first, scores[i].second) + 0.000000001 * (scores[i].first + scores[i].second); // TODO other scoring functions
        if (score > bestScore)
        {
            bestScore = score;
            bestVar = i;
        }
    }

    varScore = scores[bestVar];
    // printf("c Best variable qwer: %d\n", bestVar);
    // printf("c Score: %f %f\n", varScore.first, varScore.second);
    if (bestVar == 0)
    {
        printf("Error: no variable available\n");
        EXIT_UNWANTED_STATE
    }
    return bestVar;
}

// Extend cube trail with the given literals, propagate and return number of propagations; returns the the result of propagation
int GraphSolver::updateCubeTrail(const vector<int> lits, cubing_state_t &cubingState, int &numPropagated)
{
    LOG(LOG_LEVEL_DEBUG, "Update cube trail");
    // for (int lit : lits)
    // {
    //     assert(!isAssigned[abs(lit)]);
    //     assert(!IS_FIXED(abs(lit)));
    //     assert(active(abs(lit)));
    // }

    vector<int> assumptions;
    for (auto lit : lits)
        assumptions.push_back(lit);

    for (auto trailLits : cubingState.decisions)
    {
        // all decisions
        for (auto lit : trailLits)
            assumptions.push_back(lit);
    }

    for (auto lit : assumptions)
        CaDiCaL::Solver::assume(lit);

    clock_t start = clock();
    int res = this->propagate();
    fflush(stdout);
    LOG(LOG_LEVEL_DEBUG, "Time for propagation: " + std::to_string((double)(clock() - start) / CLOCKS_PER_SEC));
    if (res == 20) // UNSAT
    {
        // printf("Number of checks before conflict %d\n", countChecks);
        reset_assumptions();
        numPropagated = vars() + 1; // TODO just a large constant
    }
    else if (res == 10) // SAT
    {
        reset_assumptions();

        cubingState.cubeTrail.push_back(lits);
        cubingState.decisions.push_back(lits);
        for (int lit : lits)
            cubingState.isAssigned[abs(lit)] = true;

        // go over all variables and check truth value in model
        for (int i = 1; i <= this->vars(); i++)
        {
            if (!cubingState.isAssigned[i] && !IS_FIXED(i) && active(i))
                cubingState.cubeTrail.back().push_back(this->val(i));
        }

        for (auto l : cubingState.cubeTrail.back())
            cubingState.isAssigned[abs(l)] = true;

        numPropagated = (int)cubingState.cubeTrail.back().size();
    }
    else
    {
        // get number of assigned literals
        vector<int> entrailed_literals;
        get_entrailed_literals(entrailed_literals);
        reset_assumptions();

        cubingState.cubeTrail.push_back(lits);
        cubingState.decisions.push_back(lits);
        for (int lit : lits)
            cubingState.isAssigned[abs(lit)] = true;

        for (auto l : entrailed_literals)
        {
            if (!cubingState.isAssigned[abs(l)] && !IS_FIXED(abs(l)))
            {
                cubingState.cubeTrail.back().push_back(l);
                numPropagated++;
            }
        }

        for (auto l : cubingState.cubeTrail.back())
            cubingState.isAssigned[abs(l)] = true;

        numPropagated = (int)cubingState.cubeTrail.back().size();
    }
    return res;
}

void GraphSolver::undoCubeTrail(cubing_state_t &cubingState)
{
    for (auto lit : cubingState.cubeTrail.back())
    {
        if (!(IS_FIXED(abs(lit))))
        {
            cubingState.isAssigned[abs(lit)] = false;
        }
    }
    cubingState.cubeTrail.pop_back();
    cubingState.decisions.pop_back();
}

// traverse clause
class RelevantVarFinder : public CaDiCaL::ClauseIterator
{

public:
    std::vector<int> assignment;     // the current assignment of the variable
    std::vector<int> occurenceCount; // number of occurences as positive literal
    std::vector<int> occurenceCountPos;
    std::vector<int> occurenceCountNeg;

    RelevantVarFinder(int numVars, vector<vector<int>> trail)
    {
        assignment.resize(numVars + 1, 0);
        occurenceCount.resize(numVars + 1, 0);
        occurenceCountPos.resize(numVars + 1, 0);
        occurenceCountNeg.resize(numVars + 1, 0);

        for (auto lits : trail)
            for (auto lit : lits)
                assignment[abs(lit)] = lit;
    }

    bool clause(const std::vector<int> &clause)
    {
        // check if already satisfied
        for (auto lit : clause)
        {
            if (assignment[abs(lit)] == lit)
                return true; // ignore satisfied clauses
        }

        for (auto lit : clause)
        {
            occurenceCount[abs(lit)]++;
            if (lit > 0)
                occurenceCountPos[abs(lit)]++;
            else
                occurenceCountNeg[abs(lit)]++;
        }
        return true;
    };
};

bool GraphSolver::lookahead(cubing_state_t &cubingState, const bool repeatOnConflict, vector<std::pair<int, int>> &numPropagated)
{
    // this->printStats();

    LOG(LOG_LEVEL_DEBUG, "Lookahead");
    inLookaheadState = true;

    clock_t timeInSolver = 0;

    vector<int> relevantVars;
    int numVars = vars();
    if (config.lookaheadOnlyEdgeVars)
        numVars = graphHandler->getNumVariables();

    // traverse clauses to check which variables are really relevant

    // RelevantVarFinder relevantVarFinder(numVars, cubeTrail);
    // CaDiCaL::Solver::traverse_clauses(relevantVarFinder, false, false);

    // for (int i = 1; i <= numVars; i++)
    // {
    //     if (IS_FIXED(i) || relevantVarFinder.occurenceCount[i] <= 1 || !active(i) || isAssigned[i])
    //         continue;
    //     relevantVars.push_back(i);
    // }

    for (int i = 1; i <= numVars; i++)
    {
        if (IS_FIXED(i) || !active(i) || cubingState.isAssigned[i])
            continue;
        relevantVars.push_back(i);
    }

    LOG(LOG_LEVEL_DEBUG, "Number of relevant variables: " + std::to_string(relevantVars.size()));

    int loopCounter = 0;
    while (true)
    {
        int testCounter = 0;
        loopCounter++;
        LOG(LOG_LEVEL_DEBUG, "Start checking variables in lookahead: round" + std::to_string(loopCounter));
        LOG(LOG_LEVEL_DEBUG, "Number of relevant variables: " + std::to_string(relevantVars.size()));

        int countChecks = 0;
        int countConflicts = 0;

        for (int i : relevantVars)
        {

            if (IS_FIXED(i) || !active(i) || cubingState.isAssigned[i])
            {
                // just to be sure
                numPropagated[i].first = 0;
                numPropagated[i].second = 0;
                continue;
            }

            // printf("Testing variable %d\n", i);

            for (auto lit : {i, -i})
            {
                testCounter++;

                countChecks++;

                int r;
                clock_t timeInSolverStart = clock();
                if (lit == i)
                    r = updateCubeTrail({lit}, cubingState, numPropagated[i].first);
                else
                    r = updateCubeTrail({lit}, cubingState, numPropagated[i].second);

                timeInSolver += clock() - timeInSolverStart;

                if (r == 20) // UNSAT
                {
                    // printf("Conflict\n");
                    countConflicts++;
                    int numProp;
                    int rNegated = updateCubeTrail({-lit}, cubingState, numProp); // check whether negated literal is still valid

                    if (rNegated == 20)
                    {
                        inLookaheadState = false;
                        LOG(LOG_LEVEL_DEBUG, "Total time for solver: " + std::to_string((double)timeInSolver / CLOCKS_PER_SEC));
                        return false;
                    }
                    else
                    {
                        undoCubeTrail(cubingState); // undo the last assignment
                    }

                    if (repeatOnConflict)
                    {
                        // have to update the trail so there is progress and the same variable is not checked again, because officially not assigned yet.
                        auto decisions = cubingState.decisions.back();
                        undoCubeTrail(cubingState);
                        int numProp;
                        int res2 = updateCubeTrail(decisions, cubingState, numProp);
                        if (res2 == 20)
                        {
                            throw std::runtime_error("Literal already tested and shouldn't lead to a conflict");
                        }

                        numPropagated[i].first = 0;
                        numPropagated[i].second = 0;
                        break; // no point in checking other lit and also problems because already assigned
                    }

                    // printf("Conlflicting lit: %d; %d\n", lit, i);
                    // printf("Number of checks before conflict %d\n", countChecks);
                }
                else
                {
                    undoCubeTrail(cubingState);
                }
            }
        }

        // printf("Count checks %d\n", countChecks);
        // printf("Count conflicts %d\n", countConflicts);

        if (!repeatOnConflict)
            break;
        else if (countConflicts == 0)
            break;

        LOG(LOG_LEVEL_DEBUG, "Test counter" + std::to_string(testCounter));
    }
    inLookaheadState = false;
    LOG(LOG_LEVEL_DEBUG, "Total time for solver: " + std::to_string((double)timeInSolver / CLOCKS_PER_SEC));
    return true;
}

bool GraphSolver::checkCutoff(const vector<vector<int>> &cubeTrail, int assignmentCutoff)
{
    if (countAssignedOrInactive(cubeTrail) >= assignmentCutoff)
    {
        // printf("Cutoff reached\n");
        printCube(cubeTrail, fixedLiterals, config.cubesOnlyDecisions);
        return false;
    }
    return true;
}

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

void GraphSolver::generateCubesRec(cubing_state_t &cubingState, int assignmentCutoff)
{

    if (!checkCutoff(cubingState.cubeTrail, assignmentCutoff))
        return;

    // pick branching variable TODO pick more recenable branching variable
    int branchingVar = 0;

    bool simpleVariableSelection = false;
    bool useLookahead = true;
    if (simpleVariableSelection)
    {
        for (int i = 1; i < (int)cubingState.isAssigned.size(); i++)
        {
            if (!cubingState.isAssigned[i] && !(IS_FIXED(i)) && active(i)) // not assigned and not fixed
            {
                branchingVar = i;
                break;
            }
        }
    }
    else if (useLookahead)
    {
        vector<std::pair<int, int>> scores(vars() + 1, {0, 0}); // reset

        clock_t start = clock();
        // ExtendedClauseCounter counter(100000);
        // CaDiCaL::Solver::traverse_clauses(counter, true, false);
        // printf("Number of clauses before: %d\n", counter.count);
        // int before = counter.count;
        bool r = lookahead(cubingState, true, scores);
        // ExtendedClauseCounter counter2(100000);
        // CaDiCaL::Solver::traverse_clauses(counter2, true, false);
        // printf("Number of clauses after: %d\n", counter2.count);
        // if (counter2.count > before)
        // {
        //     printf("Error: More clauses after lookahead\n");
        //     EXIT_UNWANTED_STATE
        // }
        LOG(LOG_LEVEL_DEBUG, "Time for lookahead: " + std::to_string((double)(clock() - start) / CLOCKS_PER_SEC));

        if (!r)
        {
            return; // conflict under current assumptions
        }
        if (!checkCutoff(cubingState.cubeTrail, assignmentCutoff))
            return;

        std::pair<int, int> varScore;
        branchingVar = getBestVar(scores, varScore);
        if (varScore.first < varScore.second)
            branchingVar = -branchingVar; // pick the one which is more likely to quickly lead to a conflict
    }
    else
    {
        EXIT_UNWANTED_STATE
    }

    if (!branchingVar)
    {
        printf("No branching variable left");
        EXIT_UNWANTED_STATE
    }

    LOG(LOG_LEVEL_DEBUG, "Branching variable: " + std::to_string(branchingVar));
    // printf("c Best var: %d\n", branchingVar);

    for (auto branchingLit : {branchingVar, -branchingVar})
    {
        int numPropagated;
        int res = this->updateCubeTrail({branchingLit}, cubingState, numPropagated);
        if (res == 20) // UNSAT
        {
            continue;
        }
        else if (res == 10) // SAT
        {
            printCube(cubingState.cubeTrail, fixedLiterals, config.cubesOnlyDecisions); // also solution as cube
            printf("Found a solution during cubing!!!!\n");
        }

        generateCubesRec(cubingState, assignmentCutoff);
        undoCubeTrail(cubingState);
    }
}

// print as python list
void printGame(const vector<vector<int>> &cubeTrail)
{
    printf("[");
    for (int i = 0; i < (int)cubeTrail.size(); i++)
    {
        printf("[");
        for (int j = 0; j < (int)cubeTrail[i].size(); j++)
        {
            printf("%d", cubeTrail[i][j]);
            if (j < (int)cubeTrail[i].size() - 1)
                printf(", ");
        }
        printf("]");
        if (i < (int)cubeTrail.size() - 1)
            printf(", ");
    }
    printf("]\n");
}

int GraphSolver::countAssignedOrInactive(const vector<vector<int>> &cubeTrail)
{
    int count = 0;
    for (auto lits : cubeTrail)
        for (auto lit : lits)
        {
            if (active(abs(lit)))
                count++;
        }

    count += vars() - active(); // inactive variables
    return count;
}

bool GraphSolver::createGameRec(cubing_state_t &cubingState, int recDepth, int &selectedBranchingLit, int &finalScore)
{
    vector<std::pair<int, int>> scores(vars() + 1, {0, 0}); // reset;
    clock_t start = clock();
    bool r = lookahead(cubingState, false, scores);
    LOG(LOG_LEVEL_DEBUG, "Time for lookahead: " + std::to_string((double)(clock() - start) / CLOCKS_PER_SEC));
    if (!r)
    {
        return false; // conflict under current assumptions
    }

    if (std::all_of(scores.begin() + 1, scores.end(), [](const std::pair<int, int> &score)
                    { return score.first == 0 && score.second == 0; }))
    {
        return false; // all variables are assigned
    }

    std::pair<int, int> varScore;
    int branchingLit = getBestVar(scores, varScore);

    if (recDepth == 0)
    {
        if (varScore.first < varScore.second)
            branchingLit = -branchingLit; // pick the one which leads to a longer game, i.e., a longer path on the search tree
        selectedBranchingLit = branchingLit;
        finalScore = std::min(varScore.first, varScore.second);
        return true;
    }

    int LARGE_CONSTANT = 10000000;
    finalScore = LARGE_CONSTANT;

    selectedBranchingLit = branchingLit;

    for (auto lit : {branchingLit, -branchingLit})
    {
        printf("Checkout: %d %d\n", lit, recDepth - 1);

        int numPropagated;
        this->inLookaheadState = true;
        int res = this->updateCubeTrail({lit}, cubingState, numPropagated);
        this->inLookaheadState = false; // avoid adding additional clauses by propagator
        if (res != 0)                   // not an open case anymore, either 10 or 20
        {
            if (res == 20)
            {
                // undoCubeTrail(isAssigned, cubeTrail);
                continue;
            }

            if (res == 10)
            {
                undoCubeTrail(cubingState);
                continue;
            }
            printf("Result: %d\n", res);
            printf("Trail: ");
            for (auto lits : cubingState.cubeTrail)
                for (auto lit : lits)
                    printf("%d ", lit);
            printf("\n");
            EXIT_UNWANTED_STATE; // should be an open case
        }

        int selectedBranchingLitRec;
        int finalScoreRec;

        if (createGameRec(cubingState, recDepth - 1, selectedBranchingLitRec, finalScoreRec))
        {
            int summedScore = finalScoreRec + (lit > 0 ? varScore.first : varScore.second);
            if (summedScore < finalScore)
            {
                finalScore = summedScore;
                selectedBranchingLit = lit;
            }
        }
        else
        {
            // noting todo
        }

        undoCubeTrail(cubingState);
    }

    if (finalScore == LARGE_CONSTANT)
    {
        return false;
    }

    return true;
}

void GraphSolver::createGame(float randomize, vector<int> cube, int recLvl)
{
    if (randomize)
        srand(time(NULL));

    this->simplify(); // helps to drastically speed up the search
    this->set("nolearn", 1); // deactivate learning on conflicts when using propagate

    for (auto lit: cube)
    {
        add(lit);
        add(0);
    }

    // printf("Number of vars: %d\n", this->vars());
    cubing_state_t cubingState = {};
    cubingState.isAssigned = vector<bool>(this->vars() + 1, false);

    int numPropagated;
    int res = updateCubeTrail(vector<int>(), cubingState, numPropagated); // get root level assigned variables
    if (res == 20)                                               // UNSAT
    {
        printf("[]\n"); // empty trail
        printf("Warning: assumptions in conflict with encoding\n");
        return;
    }

    LOG(LOG_LEVEL_INFO, "Number of variables: " + std::to_string(vars()));
    LOG(LOG_LEVEL_INFO, "Number of active variables: " + std::to_string(active()));
    fflush(stdout);

    // printf("c Initial trail: ");
    // for (auto lit : cubeTrail.back())
    //     printf("%d ", lit);
    // printf("\n");

    while (true)
    {
        int selectedBranchingLit;
        int finalScore;
        if (!createGameRec(cubingState, recLvl, selectedBranchingLit, finalScore))
        {
            printGame(cubingState.cubeTrail);
            return; // conflict under current assumptions
        }

        // flip biased on randomize
        if (randomize > 0.0)
        {
            if (rand() % 10000 < randomize * 10000)
                selectedBranchingLit = -selectedBranchingLit;
        }

        printf("Branching lit: %d; score %d\n", selectedBranchingLit, finalScore);
        fflush(stdout);

        int numPropagated;
        int res = this->updateCubeTrail({selectedBranchingLit}, cubingState, numPropagated);
        if (res != 0) // not an open case anymore, either 10 or 20
        {
            printGame(cubingState.cubeTrail);
            return;
        }

        // also abbort if all edge variables are assigned or fixed or inactive
        int numEdgeVars = graphHandler->getNumVariables();
        bool allEdgeVarsAssignedFixedOrInactive = true;
        for (int i = 1; i <= numEdgeVars; i++)
        {
            if (!cubingState.isAssigned[i] && !IS_FIXED(i) && active(i))
            {
                allEdgeVarsAssignedFixedOrInactive = false;
                break;
            }
        }
        if (allEdgeVarsAssignedFixedOrInactive)
        {
            printGame(cubingState.cubeTrail);
            return;
        }
    }

    return;
}