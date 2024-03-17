#ifndef CADICAL_SOLVER_INTERFACE_H
#define CADICAL_SOLVER_INTERFACE_H

#include "sms.hpp"
#include "cadical.hpp"
#include <deque>

using std::deque;

class CadicalSolver : public GraphSolver, public CaDiCaL::ExternalPropagator
{
private:
    bool redundant;

    bool changeInTrail = true; // checks whether the trail has changed since the last propagation step

    int highestEdgeVariable;

    vector<vector<int>> clauses;
    int incrementalMode = false; // if true solver has finished and clauses are added by the normal "incremental interface", i.e., adding clauses without observed variables is possible

    deque<vector<int>> current_trail; // for each decision lvl store the assigned literals (only positive version)
    vector<bool> isFixed;             // isFixed[v] is true if the truth value of this variable is fixed

    vector<vector<int>> literal2clausePos; // for each edge variable store clause which was used the last time.
    vector<vector<int>> literal2clauseNeg; // for each negation of an edge variable

    void init(SolverConfig config, cnf_t &cnf);

public:
    CaDiCaL::Solver *solver;
    CadicalSolver(SolverConfig config);
    CadicalSolver(SolverConfig config, cnf_t &cnf);
    ~CadicalSolver()
    {
        solver->disconnect_external_propagator();
        delete solver;
    }

    bool solve(vector<int> assumptions);
    bool solve(vector<int> assumptions, int timeout);
    void printFullModel(void);

    /* API function to return and internally block the next solution
     * the returned data has the following format:
     *    
     *    // g = getNextGraph()
     *    *g is an int that holds the number of edges, m
     *    g[1]-g[2] ... g[2m-1]-g[2m] are the edges of the graph
     *
     * the data is stored in last_graph, and the returned pointer
     * points to its beginning
     *   
     */
    int* getNextGraph(vector<int> assumptions);
    vector<int> last_graph;

public:
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
        currentAssignment[absLit] = lit > 0 ? truth_value_true : truth_value_false;
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
                    currentAssignment[l] = truth_value_unknown;
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
        if (config.checkSolutionInProp)
        {
            return check();
        }
        return true;
    }

    bool check_solution()
    {
        // printf("Check whether there are fixed variables\n");
        // int nFixed = 0;
        // for (int i = 1; i < nextFreeVariable; i++)
        // {
        //     if (solver->fixed(i) !=  0)
        //     {
        //         // printf("Variable %d is fixed\n", i);
        //         nFixed++;
        //     }
        // }
        // printf("Number of fixed variables: %d\n", nFixed);
        if (!config.checkSolutionInProp)
        {
            incrementalMode = true;
            vector<int> currentModel;
            this->model = &currentModel; // TODO have to extract current model because the one from cb_check is deleted
            for (int i = 1; i < nextFreeVariable; i++)
            {
                currentModel.push_back(solver->val(i));
            }
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
            if (currentAssignment[absLit] == truth_value_unknown)
            {
                nUnknown++;
                unassigned = l;
            }
            else if (currentAssignment[absLit] == truth_value_true && l > 0)
                return 0; // already satisfied
            else if (currentAssignment[absLit] == truth_value_false && l < 0)
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

extern "C" {
  int* next_solution(void* sms_solver);
  void* create_solver(int vertices);
  void destroy_solver(void* sms_solver);
  void add_literal(void* sms_solver, int lit);
}

#endif
