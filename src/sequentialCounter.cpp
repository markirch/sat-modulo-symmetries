#include "sequentialCounter.hpp"

inline void addClause(CaDiCaL::Solver *s, clause_t c)
{
    for (auto l : c)
        s->add(l);
    s->add(0);
}

// -1 is interpreted as don't care. Returns the variables that indicate the number of true values in vars
vector<int> sequentialCounter(CaDiCaL::Solver *solver, const vector<int> &vars, const int countUpTo, int &varCount, const int atLeast, const int atMost)
{
    int n = vars.size();
    vector<vector<int>> counterVariables(n, vector<int>(countUpTo));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < countUpTo; j++)
            counterVariables[i][j] = ++varCount; // Create new variables

    counterVariables[0][0] = vars[0];
    for (int i = 1; i < countUpTo; i++)
        add_clause(solver, {-counterVariables[0][i]}); // At most one at the beginning

    for (int i = 0; i < n - 1; i++)
    {
        add_clause(solver, {-vars[i + 1], counterVariables[i + 1][0]}); // If there is an element then there is at least one element
        for (int j = 0; j < countUpTo; j++)
        {
            add_clause(solver, {-counterVariables[i][j], counterVariables[i + 1][j]});              // At least as many
            add_clause(solver, {counterVariables[i][j], vars[i + 1], -counterVariables[i + 1][j]}); // The same if element is not present

            if (j < countUpTo - 1)
            {
                add_clause(solver, {-counterVariables[i][j], -vars[i + 1], counterVariables[i + 1][j + 1]}); // One more element
                add_clause(solver, {counterVariables[i][j], -counterVariables[i + 1][j + 1]});               // At most one more
            }
        }
    }

    if (atMost != -1)
        for (int i = 0; i < n - 1; i++)
            add_clause(solver, {-counterVariables[i][atMost - 1], -vars[i + 1]}); // If maximum reached, no more true variables
    if (atLeast != -1)
        add_clause(solver, {counterVariables[n - 1][atLeast - 1]});

    vector<int> result(countUpTo);
    for (int j = 0; j < countUpTo; j++)
        result[j] = counterVariables[n - 1][j];
    varCount++; // little wasteful but potentially avoids some errors if function is not used correctly
    return result;
}