#ifndef CLINGO_SOLVER_INTERFACE_H
#define CLINGO_SOLVER_INTERFACE_H

#include <clingo.h>
#include "sms.hpp"

class ClingoSolver : public GraphSolver
{
public:
    ClingoSolver(SolverConfig &config, cnf_t &cnf);
    ~ClingoSolver() {}
    clingo_propagate_control_t *propagate_control;
    bool clauseAddable;      // true if clause can be added in the current propagation step; false if a conflicting clause was already added
    bool addStoredClauses(); // add as many clauses stored as possible; return true if no conflicting clause was added

    vector<lit_t> variablesOfSubgraphMapped; // just a list of variables representing the subgraph
protected:                                   // virtual classes from common interface
    bool solve(vector<int> assumptions);
    bool solve(vector<int> assumptions, int timeout);
    adjacency_matrix_t getAdjacencyMatrix();
    void addClause(const vector<lit_t> &clause, bool redundant);

private:
    clingo_control_t *ctl;
    vector<clingo_atom_t> atoms;                  // for each integer in SAT encoding corresponding atom, so atom[i] gives the atom representing the variable i
    vector<lit_t> variables2solverVariables;      // mapping from variables from the original encoding to clingo literals during solving. Only for the once mapped in initialization
    vector<vector<clingo_literal_t>> mappedEdges; //
    vector<vector<clingo_literal_t>> mappedEdgesIntersectionGraph; //
    vector<vector<vector<clingo_literal_t>>> mappedTriangles;

    vector<vector<lit_t>> redundandentClauses;
    vector<vector<lit_t>> irredundentClauses;

public:
    // the functions for the propagator (class functions are not possible, since the cann't be given to clingo as argument)
    static bool propagate_clingo(clingo_propagate_control_t *control, const clingo_literal_t *, size_t, ClingoSolver *s);
    static bool undo_clingo(clingo_propagate_control_t *control, const clingo_literal_t *, size_t, ClingoSolver *s);
    static bool check_clingo(clingo_propagate_control_t *control, ClingoSolver *s);
    static bool init_clingo(clingo_propagate_init_t *init, ClingoSolver *s);

public:
    void init(clingo_propagate_init_t *ctl);
};

#endif
