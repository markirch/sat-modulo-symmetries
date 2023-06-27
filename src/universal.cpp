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