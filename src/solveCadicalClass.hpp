#ifndef CADICAL_SOLVER_INTERFACE_H
#define CADICAL_SOLVER_INTERFACE_H

#include "solver.hpp"
#include "cadical.hpp"
#include <bits/stdc++.h>
#include <algorithm>

class CadicalSolver : public GraphSolver, public CaDiCaL::ExternalPropagator
{
private:
    bool redundant;

    bool changeInTrail = true; // checks whether the trail has changed since the last propagation step

    int highestEdgeVariable;

    vector<vector<int>> clauses;
    int highestVariable;
    int incrementalMode = false; // if true solver has finished and clauses are added by the normal "incremental interface", i.e., adding clauses without observed variables is possible

    deque<vector<int>> current_trail;       // for each decision lvl store the assigned literals (only positive version)
    vector<truth_value_t> currentAssigment; // currentAssigment[v] gives the truthvalue of variable v (if observed)
    vector<bool> isFixed;                   // isFixed[v] is true if the truth value of this variable is fixed

    vector<vector<int>> literal2clausePos; // for each edge variable store clause which was used the last time.
    vector<vector<int>> literal2clauseNeg; // for each negation of an edge variable

    void init(configSolver config, cnf_t &cnf);

public:
    CaDiCaL::Solver *solver;
    CaDiCaL::Solver *universalSolver;
    CadicalSolver(configSolver config);
    CadicalSolver(configSolver config, cnf_t &cnf);
    ~CadicalSolver() { solver->disconnect_external_propagator(); }

    bool solve(vector<int> assumptions);
    bool solve(vector<int> assumptions, int timeout);

protected: // virtual classes from common interface
    adjacency_matrix_t getAdjacencyMatrix()
    {
        // printf("Trail: ");
        // for (auto lit: *current_trail)
        //     printf("%d ", lit);
        // printf("\n");
        if (!incrementalMode)
        {
            adjacency_matrix_t matrix(vertices, vector<truth_value_t>(vertices, truth_value_unknown));
#ifndef DIRECTED
            for (int i = 0; i < vertices; i++)
                for (int j = i + 1; j < vertices; j++)
                    matrix[i][j] = matrix[j][i] = currentAssigment[edges[i][j]];
#else
            for (int i = 0; i < vertices; i++)
                for (int j = 0; j < vertices; j++)
                    if (i != j)
                        matrix[i][j] = currentAssigment[edges[i][j]];

#endif
            // printFullMatrix = true;
            // printAdjacencyMatrix(matrix);
            return matrix;
        }
        else
        {
            adjacency_matrix_t matrix(vertices, vector<truth_value_t>(vertices, truth_value_unknown));
#ifndef DIRECTED
            for (int i = 0; i < vertices; i++)
                for (int j = i + 1; j < vertices; j++)
                    matrix[i][j] = matrix[j][i] = solver->val(edges[i][j]) > 0 ? truth_value_true : truth_value_false;
#else
            for (int i = 0; i < vertices; i++)
                for (int j = 0; j < vertices; j++)
                    if (i != j)
                        matrix[i][j] = solver->val(edges[i][j]) > 0 ? truth_value_true : truth_value_false;

#endif
            // printFullMatrix = true;
            // printAdjacencyMatrix(matrix);
            return matrix;
        }
    }

    vector<adjacency_matrix_t> getAdjacencyMatrixMultiple()
    {
        // printf("Trail: ");
        // for (auto lit: *current_trail)
        //     printf("%d ", lit);
        // printf("\n");
        int nMatrices = config.edgesMultiple.size();
        vector<adjacency_matrix_t> matrices(nMatrices, adjacency_matrix_t(vertices, vector<truth_value_t>(vertices, truth_value_unknown)));
#ifndef DIRECTED
        for (int n = 0; n < nMatrices; n++)
            for (int i = 0; i < vertices; i++)
                for (int j = i + 1; j < vertices; j++)
                    matrices[n][i][j] = matrices[n][j][i] = currentAssigment[config.edgesMultiple[n][i][j]];
#else
        EXIT_UNWANTED_STATE // not supported yet
#endif
        // printFullMatrix = true;
        // printAdjacencyMatrix(matrix);
        return matrices;
    }

    vector<truth_value_t> &getCurrentAssignemnt() { return currentAssigment; }

    vector<vector<truth_value_t>> getStaticPartition()
    {
        vector<vector<truth_value_t>> partition(vertices, vector<truth_value_t>(vertices, truth_value_unknown));
        for (int i = 0; i < vertices; i++)
            for (int j = i + 1; j < vertices; j++)
            {
                partition[i][j] = partition[j][i] = currentAssigment[config.staticPartition[i][j]];
            }
        return partition;
    }

    void addClause(const vector<lit_t> &clause, bool)
    {

        // if (sym_breaking_clause.size() != 0)
        //    return; // EXIT_UNWANTED_STATE
        // printf("Number of literals: %ld, Add the following clause:", clause.size());
        // for (auto lit : clause)
        //     printf("%d ", lit);
        // printf("\n");
        if (!incrementalMode)
            this->clauses.push_back(clause);
        else
        {
            if (config.addedClauses)
            {
                for (auto l : clause)
                    fprintf(config.addedClauses, "%d ", l);
                fprintf(config.addedClauses, "0\n");
            }
            // use incremental interface
            for (auto l : clause)
                solver->add(l);
            solver->add(0);
        }
    }

public:
    void notify_assignment(int lit, bool is_fixed)
    {
        changeInTrail = true;
        int absLit = abs(lit);
        currentAssigment[absLit] = lit > 0 ? truth_value_true : truth_value_false;
        this->isFixed[absLit] = is_fixed;
        if (!is_fixed)
            current_trail.back().push_back(absLit);
    }

    void notify_new_decision_level()
    {
        current_trail.push_back(vector<int>());
    }

    void notify_backtrack(size_t new_level)
    {
        while (current_trail.size() > new_level + 1)
        {
            auto last = current_trail.back();
            for (int l : last)
            {
                if (!isFixed[l])
                    currentAssigment[l] = truth_value_unknown;
            }
            current_trail.pop_back();
        }
    }

    // currently not checked in propagator but with the normal incremental interface to allow adding other literals or even new once.
    bool cb_check_found_model(const std::vector<int> &model)
    {
        this->model = &model;
        if (!clauses.empty())
            return false; // EXIT_UNWANTED_STATE only do check if there isn't another clause to add before
        // this->current_trail = &model;
        if (config.chechSolutionInProp)
        {
            return check();
        }
        return true;
    }

    bool check_solution()
    {
        if (!config.chechSolutionInProp)
        {
            incrementalMode = true;
            bool res = check();
            incrementalMode = false;
            return res;
        }
        return true;
    }

    bool cb_has_external_clause()
    {
        // PRINT_CURRENT_LINE
        // if no clause, then check whether a clause could be added. If already a clause present then just return clause.
        // if propagation is done in other function then not compute clauses here
        if (clauses.empty() && changeInTrail && !config.propagateLiteralsCadical)
        {
            changeInTrail = false;
            propagate();
        }

        // printf("Check for external clause: %ld\n", sym_breaking_clause.size());
        return !clauses.empty();
    }

    int cb_add_external_clause_lit()
    {
        // PRINT_CURRENT_LINE
        // printf("Call: Add external clause\n");
        vector<int> &lastClause = clauses[clauses.size() - 1];
        if (lastClause.empty())
        {
            clauses.pop_back(); // delete last clause
            // printf(" end clause\n");
            if (config.addedClauses)
                fprintf(config.addedClauses, "0\n");
            return 0;
        }
        else
        {
            // printf("Add external clause\n");
            int lit = lastClause.back();
            lastClause.pop_back();
            // printf("%d ", lit);
            if (config.addedClauses)
                fprintf(config.addedClauses, "%d ", lit);
            return lit;
        }
    }

    // functions need to be defined
    int cb_decide() { return 0; }

    int cb_propagate()
    {
        if (!config.propagateLiteralsCadical)
            return 0;

        if (!changeInTrail)
            return 0;

        changeInTrail = false;
        propagate();

        if (clauses.empty())
            return 0;

        auto lastClause = clauses.back();
        assert(!lastClause.empty());
        // find unassigned literal otherwise take last one; first check if clause is unit
        int nUnknown = 0;
        int unassigned = 0;
        for (auto l : lastClause)
        {
            auto absLit = abs(l);
            if (currentAssigment[absLit] == truth_value_unknown)
            {
                nUnknown++;
                unassigned = l;
            }
            else if (currentAssigment[absLit] == truth_value_true && l > 0)
                return 0; // already satisfied
            else if (currentAssigment[absLit] == truth_value_false && l < 0)
                return 0; // already satisfied
        }

        if (nUnknown == 1)
        {
            // notify_assignment(unassigned, false); // push back the clause to the current trail
            clauses.pop_back(); // delete last clause
            if (unassigned > 0)
                literal2clausePos[abs(unassigned)] = lastClause;
            else
                literal2clauseNeg[abs(unassigned)] = lastClause;
            // PRINT_CURRENT_LINE
            return unassigned;
        }
        return 0;
    }

    int cb_add_reason_clause_lit(int plit)
    {
        // PRINT_CURRENT_LINE
        if (plit > 0)
        {
            if (literal2clausePos[abs(plit)].empty())
            {
                if (config.addedClauses)
                    fprintf(config.addedClauses, "0\n");
                return 0;
            }
            auto l = literal2clausePos[abs(plit)].back();
            literal2clausePos[abs(plit)].pop_back();
            if (config.addedClauses)
                fprintf(config.addedClauses, "%d ", l);
            return l;
        }
        else
        {
            if (literal2clauseNeg[abs(plit)].empty())
            {
                if (config.addedClauses)
                    fprintf(config.addedClauses, "0\n");
                return 0;
            }
            auto l = literal2clauseNeg[abs(plit)].back();
            literal2clauseNeg[abs(plit)].pop_back();
            if (config.addedClauses)
                fprintf(config.addedClauses, "%d ", l);
            return l;
        }
    };
};

#endif
