#include "universal.hpp"

void UniversalChecker::checkProperty(const adjacency_matrix_t &, const vector<int> &model, int &nextFreeVariable)
{
    for (auto var : forAllAsumptions)
        if (model[var - 1] > 0)
            universalSolver->assume(var);
        else
            universalSolver->assume(-var);
    int res = universalSolver->solve();
    if (res == 10)
    {
        clause_t clauseForAll;
        vector<clause_t> clausesForSolver;
        for (auto clause : forAllCNF)
        {
            vector<int> clauseReplaced; // clause where we replace assigned
            bool satisfied = false;     // if clause already true
            for (auto lit : clause)
            {
                if (std::find(forAllAsumptions.begin(), forAllAsumptions.end(), abs(lit)) != forAllAsumptions.end())
                {
                    // element is an assumption
                    clauseReplaced.push_back(lit);
                }
                else
                {
                    // literal is universal
                    if (universalSolver->val(abs(lit)) == lit)
                    {
                        satisfied = true;
                        break;
                    }
                }
            }

            if (!satisfied)
            {
                // construct clause; i.e., at least one clause must not be satisfied for the universal susbstitution
                if (clauseReplaced.size() == 1)
                {
                    clauseForAll.push_back(-clauseReplaced[0]);
                }
                else
                {
                    // new variable
                    int newVar = nextFreeVariable++; // if the new variable is true, then all literals in the clause are false
                    for (auto lit : clauseReplaced)
                    {
                        clausesForSolver.push_back({-newVar, -lit});
                        clauseForAll.push_back(newVar);
                    }
                }
            }
        }
        clausesForSolver.push_back(clauseForAll); // this clause is the one which says that one clause must be false in the other one
        throw clausesForSolver;
    }
}

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>

void UniversalCheckerQCIR::checkProperty(const adjacency_matrix_t &, const vector<int> &model, int &nextFreeVariable)
{
    // only applicable if directly translated from QCIR with Teitsin (without polarity optimizations)
    // clock_t startTotal = clock();
    int outputGate = gateVariables[indexOutputGate]; // TODO ensure that first  (or last and change) getvariable is output gate

    // printf("Modelsize: %ld\n", model.size());
    // printf("Existential model: ");
    for (auto var : forAllAsumptions)
    {
        if ((int)model.size() <= var - 1)
        {
            printf("Error: variable assignment is not in scope\n");
            EXIT_UNWANTED_STATE
        }
        if (model[var - 1] > 0) // TODO if not present just assume true
            universalSolver->assume(var);
        else
            universalSolver->assume(-var);

        // printf("%d ", model[var - 1]);
    }
    // printf("\n");

    // clock_t start = clock();
    int res = universalSolver->solve();
    // print time spend in last function
    // clock_t end = clock();
    // double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
    // printf("Time: %f\n", elapsed_secs);
    if (res == 10)
    {
        bool learnFunctionInCurrentStep = false;
        vector<int> universalVariablesNew;
        vector<clause_t> clausesForSolver;

        // if (functionLearning)
        // {
        //     vector<bool> existentialAssignment;
        //     vector<bool> universalAssignment;
        //     for (auto v : forAllAsumptions)
        //         existentialAssignment.push_back(model[v - 1] > 0);
        //     for (auto v : universalVariables)
        //         universalAssignment.push_back(universalSolver->val(v) > 0);

        //     samples.push_back({existentialAssignment, universalAssignment});

        //     if (samples.size() == 40)
        //     {
        //         learnFunctionInCurrentStep = true;
        //         printf("---------------analyize---------------\n");
        //         // analyize

        //         for (int i = 0; i < (int)universalVariables.size(); i++)
        //         {
        //             vector<sample_t> treeSamples;
        //             for (auto s : samples)
        //             {
        //                 treeSamples.push_back({s.first, s.second[i]});
        //             }

        //             auto root = buildDecisionTree(treeSamples);
        //             if (root->is_leaf)
        //             {
        //                 universalVariablesNew.push_back(0); // no new variable must be assigned; truth assignment also corresponds with the last one
        //             }
        //             else
        //             {
        //                 // std::cout << "Decision tree for variable" << universalVariables[i] << std::endl;

        //                 auto v = nextFreeVariable++;
        //                 universalVariablesNew.push_back(v);
        //                 // ensure that variable has a certain truth value based on the decision tree

        //                 // printDecisionTree(root);

        //                 vector<pair<int, bool>> path;
        //                 vector<pair<vector<pair<int, bool>>, bool>> paths;
        //                 findAllPaths(root, path, paths);

        //                 for (auto p : paths)
        //                 {
        //                     clause_t clause;

        //                     for (auto f : p.first)
        //                     {
        //                         // printf("(%d,%d) ", f.first, f.second ? 1 : 0);
        //                         if (f.second)
        //                             clause.push_back(-forAllAsumptions[f.first]); // negate the path
        //                         else
        //                             clause.push_back(forAllAsumptions[f.first]);
        //                     }

        //                     if (p.second) // gives truth assignment of universal variable on this path
        //                         clause.push_back(v);
        //                     else
        //                         clause.push_back(-v);

        //                     clausesForSolver.push_back(clause);

        //                     // printf(" has value %d\n", p.second ? 1 : 0);
        //                 }
        //             }

        //             deleteDecisionTree(root);
        //         }

        //         samples.clear();
        //     }
        // }

        vector<int> gateVariablesNew;
        vector<bool> hashedInCurrentRound(gateVariables.size(), false);
        if (hashGates)
        {
            for (int i = 0; i < (int)gateVariables.size(); i++)
            {

                if (i == 0) // TODO special case, take care
                {
                    gateVariablesNew.push_back(nextFreeVariable++);
                    continue;
                }

                auto g = gateVariables[i];
                auto dep = gate2dependencies[g];
                vector<bool> state;
                for (auto v : dep)
                {
                    // printf("%d\n", v);
                    // fflush(stdout);
                    state.push_back(universalSolver->val(v) > 0);
                }
                auto m = gate2int[g];

                bool dependingOnFunction = false;
                if (m.find(state) != m.end())
                {
                    hashedInCurrentRound[i] = true;
                    gateVariablesNew.push_back(m[state]);

                    // if (learnFunctionInCurrentStep)
                    // {
                    //     for (auto v : dep)
                    //     {
                    //         auto it = std::find(universalVariables.begin(), universalVariables.end(), v);
                    //         int index = std::distance(universalVariables.begin(), it);
                    //         if (universalVariablesNew[index] != 0) // is a function
                    //         {
                    //             dependingOnFunction = true;
                    //             break;
                    //         }
                    //     }
                    // }

                    // if (!dependingOnFunction)
                    // {
                    //     // already computed at some point
                    //     // printf("Found a hashed value for gate %d\n", g);
                    //     hashedInCurrentRound[i] = true;
                    //     gateVariablesNew.push_back(m[state]);
                    // }
                }

                if (!hashedInCurrentRound[i])
                {
                    int v = nextFreeVariable++;
                    gateVariablesNew.push_back(v);
                    if (!dependingOnFunction)
                        gate2int[g][state] = v; // add new entry so it doesn't have to be calculated again
                    // printf("%d %ld, %s\n", g, state.size(), stateAsString.c_str());
                }
            }
            // end = clock();
            // double elapsed_secs = double(end - startTotal) / CLOCKS_PER_SEC;
            // printf("TimeAfter hashing: %f\n", elapsed_secs);
        }
        else
        {
            for (int i = 0; i < (int)gateVariables.size(); i++)
                gateVariablesNew.push_back(nextFreeVariable++);
        }

        for (auto clause : forAllCNF)
        {
            if (clause.size() == 1 && clause[0] == outputGate) // skip that output gate must be true
                continue;

            vector<int> clauseReplaced; // clause where we replace assigned
            bool satisfied = false;     // if clause already true

            if (hashGates)
            {
                auto gate = abs(clause[0]);
                int index = gate2index[gate];
                if (hashedInCurrentRound[index])
                {
                    // printf("Hashing had an impact\n");
                    continue;
                }
            }
            for (auto lit : clause)
            {
                if (isExistential[abs(lit)])
                {
                    // element is an assumption
                    clauseReplaced.push_back(lit);
                }
                else if (isGate[abs(lit)])
                {
                    int index = gate2index[abs(lit)];
                    if (lit > 0)
                        clauseReplaced.push_back(gateVariablesNew[index]);
                    else
                        clauseReplaced.push_back(-gateVariablesNew[index]);
                }
                else
                {
                    if (!learnFunctionInCurrentStep)
                    {
                        // literal is universal
                        if (universalSolver->val(abs(lit)) == lit)
                        {
                            satisfied = true;
                            break;
                        }
                    }
                    // else
                    // {
                    //     // replace by function
                    //     auto it = std::find(universalVariables.begin(), universalVariables.end(), abs(lit));
                    //     if (it == universalVariables.end())
                    //         EXIT_UNWANTED_STATE

                    //     int index = std::distance(universalVariables.begin(), it);
                    //     if (universalVariablesNew[index] == 0)
                    //     {
                    //         // not a function so normal case
                    //         if (universalSolver->val(abs(lit)) == lit)
                    //         {
                    //             satisfied = true;
                    //             break;
                    //         }
                    //     }
                    //     else
                    //     {
                    //         clauseReplaced.push_back(lit > 0 ? universalVariablesNew[index] : -universalVariablesNew[index]);
                    //     }
                    // }
                }
            }

            if (!satisfied)
            {
                clausesForSolver.push_back(clauseReplaced);
            }
        }

        // printf("Add clauses\n");
        clausesForSolver.push_back({-gateVariablesNew[indexOutputGate]}); //  output must be negative

        // bool checkIfConflicting = false;
        // if (checkIfConflicting)
        // {
        //     std::unordered_map<int, bool> assignment;
        //     for (auto v : forAllAsumptions)
        //         assignment[v] = model[v - 1] > 0;

        //     bool res = unitPropagation(clausesForSolver, assignment); // check if current assignment with BCP leads to a conflict, i.e., this really excludes a solution
        //     if (res)
        //     {
        //         printf("ERROR: BCP starting with assignment of the existential variables didn't lead to a conflict\n");
        //     }
        //     else
        //     {
        //         // printf("BCP lead to a conflict\n");
        //     }
        //}
        // clock_t end = clock();
        // double elapsed_secs = double(end - startTotal) / CLOCKS_PER_SEC;
        // printf("Timetotal: %f\n", elapsed_secs);
        throw clausesForSolver;
    }
    // end = clock();
    // elapsed_secs = double(end - startTotal) / CLOCKS_PER_SEC;
    // printf("Timetotal: %f\n", elapsed_secs);
}

// Function to check if a clause is satisfied
bool isClauseSatisfied(vector<int> &clause, std::unordered_map<int, bool> &assignment)
{
    for (int literal : clause)
    {
        int var = abs(literal);
        if ((literal > 0 && assignment[var]) || (literal < 0 && !assignment[var]))
        {
            return true; // The clause is satisfied
        }
    }
    return false; // The clause is not satisfied
}

// Function to perform unit propagation; return false if this leads to a conflict
bool unitPropagation(vector<vector<int>> &clauses, std::unordered_map<int, bool> &assignment)
{
    bool modified = true;
    while (modified)
    {
        modified = false;
        for (vector<int> &clause : clauses)
        {
            int unassignedCount = 0;
            int unassignedLiteral = 0;
            bool isSatisfied = false;
            for (int literal : clause)
            {
                int var = abs(literal);
                if (assignment.find(var) == assignment.end())
                {
                    unassignedCount++;
                    unassignedLiteral = literal;

                    if ((literal > 0) == assignment[var])
                    {
                        isSatisfied = true;
                        break;
                    }
                }
            }
            if (isSatisfied)
                continue;

            if (unassignedCount == 0)
            {
                return false; // Conflict detected
            }

            if (unassignedCount == 1)
            {
                int var = abs(unassignedLiteral);
                assignment[var] = (unassignedLiteral > 0);
                modified = true;
            }
        }
    }

    return true; // No conflict detected
}