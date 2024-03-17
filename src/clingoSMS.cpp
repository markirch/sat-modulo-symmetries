#include "clingoSMS.hpp"

ClingoSolver::ClingoSolver(SolverConfig config, cnf_t &cnf) : GraphSolver(config)
{
    int highestObservedVariables = *max_element(config.observedVars.begin(), config.observedVars.end());
    currentAssignment = vector<truth_value_t>(highestObservedVariables + 1, truth_value_unknown); // must be created before adding clauses because already notification can happen

    clingo_backend_t *backend;

    // create a control object and pass command line arguments
    // printf("Options %d\n", nOptionsClingo);
    // for (int i = 0; i < nOptionsClingo; i++)
    //     printf("%s\n", optionsClingo[i]);
    if (!clingo_control_new(NULL, 0, NULL, NULL, 20, &ctl) != 0)
        EXIT_UNWANTED_STATE

    // get the backend
    if (!clingo_control_backend(ctl, &backend))
        EXIT_UNWANTED_STATE
    // prepare the backend for adding rules
    if (!clingo_backend_begin(backend))
        EXIT_UNWANTED_STATE
    // get atoms for each variable.

    atoms.push_back(0); // 0 is not an atom
    for (int i = 1; i < nextFreeVariable; i++)
    {
        clingo_atom_t a;
        if (!clingo_backend_add_atom(backend, NULL, &a))
            EXIT_UNWANTED_STATE

        atoms.push_back(a);
        assert(a == atoms[i]);

        // add disjunction that either true or false
        clingo_atom_t head[] = {atoms[i]};
        if (!clingo_backend_rule(backend, true, head, 1, NULL, 0)) // choice rule, unfortunately not possible over strongly negated atoms are default negation
            EXIT_UNWANTED_STATE
    }

    variables2solverVariables = vector(nextFreeVariable, 0);

    // each clause as negated version in body (default negation used, but no problem)
    for (auto clause : cnf)
    {
        clingo_literal_t body[clause.size()];
        // map literals
        for (int i = 0; i < (int)clause.size(); i++)
            body[i] = clause[i] > 0 ? -(clingo_literal_t)atoms[clause[i]] : atoms[-clause[i]]; // not that signs are reversed, so if all true then we have a problem.

        if (!clingo_backend_rule(backend, false, NULL, 0, body, clause.size())) //  choice rule unfortunately not possible over strongly negated atoms are default negation
            EXIT_UNWANTED_STATE
    }

    // finalize the backend
    if (!clingo_backend_end(backend))
        EXIT_UNWANTED_STATE

    // register propagator
    clingo_propagator_t prop = {
        (clingo_propagator_init_callback_t)init_clingo,
        (clingo_propagator_propagate_callback_t)propagate_clingo,
        (clingo_propagator_undo_callback_t)undo_clingo, // undo_clingo
        (clingo_propagator_check_callback_t)check_clingo,
        NULL};

    if (!clingo_control_register_propagator(ctl, &prop, this, false))
        EXIT_UNWANTED_STATE

}

bool ClingoSolver::solve(vector<int> assumptions)
{
    clingo_solve_handle_t *handle;
    clingo_model_t const *model;

    vector<clingo_literal_t> mappedAssumptions; // map assumptions to atoms
    for (auto x : assumptions)
    {
        if (x > 0)
            mappedAssumptions.push_back(atoms[x]);
        else
            mappedAssumptions.push_back(-atoms[-x]);
    }

    do
    {
        if (!clingo_control_solve(ctl, clingo_solve_mode_yield, &mappedAssumptions[0], mappedAssumptions.size(), NULL, NULL, &handle))
            EXIT_UNWANTED_STATE
        if (!clingo_solve_handle_resume(handle))
            EXIT_UNWANTED_STATE
        if (!clingo_solve_handle_model(handle, &model))
            EXIT_UNWANTED_STATE

        if (check_solution()) // true if no clause was added
            return true;
    } while (model);
    return true;
}

bool ClingoSolver::solve(vector<int>, int)
{
    printf("not implemented yet\n");
    EXIT_UNWANTED_STATE
    // clingo_solve_handle_t *handle;
    // clingo_model_t const *model;

    // vector<clingo_literal_t> mappedAssumptions; // map assumptions to atoms
    // for (auto x : assumptions)
    // {
    //     if (x > 0)
    //         mappedAssumptions.push_back(atoms[x]);
    //     else
    //         mappedAssumptions.push_back(-atoms[-x]);
    // }

    // if (!clingo_control_solve(ctl, clingo_solve_mode_yield | clingo_solve_mode_async, &mappedAssumptions[0], mappedAssumptions.size(), NULL, NULL, &handle))
    //     EXIT_UNWANTED_STATE

    // if (!clingo_solve_handle_resume(handle))
    //     EXIT_UNWANTED_STATE
    // bool solvedSuccessfully;
    // clingo_solve_handle_wait(handle, timeout, &solvedSuccessfully);
    // return solvedSuccessfully;
}

bool ClingoSolver::addStoredClauses()
{
    while (redundandentClauses.size() != 0)
    {
        auto clause = redundandentClauses[redundandentClauses.size() - 1];
        bool res;
        if (!clingo_propagate_control_add_clause(propagate_control, &clause[0], clause.size(), clingo_clause_type_learnt, &res))
            EXIT_UNWANTED_STATE
        redundandentClauses.pop_back();
        if (!res)
            return false;
    }

    while (irredundentClauses.size() != 0)
    {
        auto clause = irredundentClauses[irredundentClauses.size() - 1];
        bool res;
        if (!clingo_propagate_control_add_clause(propagate_control, &clause[0], clause.size(), clingo_clause_type_static, &res))
            EXIT_UNWANTED_STATE
        irredundentClauses.pop_back();
        if (!res)
            return false;
    }
    clauseAddable = true;
    return true;
}

void ClingoSolver::addClause(const vector<lit_t> &clause, bool redundant)
{
    if (!incrementalMode)
    {
        vector<clingo_literal_t> mappedClause;
        for (auto l : clause)
        {
            if (l > 0)
                mappedClause.push_back(variables2solverVariables[l]);
            else
                mappedClause.push_back(-variables2solverVariables[-l]);
        }

        if (clauseAddable)
        {
            bool res;
            if (!clingo_propagate_control_add_clause(propagate_control, &mappedClause[0], mappedClause.size(), redundant ? clingo_clause_type_learnt : clingo_clause_type_static, &res))
                EXIT_UNWANTED_STATE
            if (!res)
                clauseAddable = false;
        }
        else
        {
            if (redundant)
                redundandentClauses.push_back(mappedClause);
            else
                irredundentClauses.push_back(mappedClause);
        }
    }
    else
    {
        clingo_backend_t *backend;
        // get the backend
        if (!clingo_control_backend(ctl, &backend))
            EXIT_UNWANTED_STATE
        // prepare the backend for adding rules
        if (!clingo_backend_begin(backend))
        {
            char const *error_message;
            if (!(error_message = clingo_error_message()))
            {
                error_message = "error";
            }
            printf("%s\n", error_message);
            // ret = clingo_error_code();
            EXIT_UNWANTED_STATE
        }

        auto maxAbsolute = abs(
            *std::max_element(clause.begin(), clause.end(),
                              [](int a, int b)
                              { return std::abs(a) < std::abs(b); }));

        for (int i = nextFreeVariable; i <= maxAbsolute; i++)
        {
            clingo_atom_t a;
            if (!clingo_backend_add_atom(backend, NULL, &a))
                EXIT_UNWANTED_STATE

            atoms.push_back(a);
            assert(a == atoms[i]);

            // add disjunction that either true or false
            clingo_atom_t head[] = {atoms[i]};
            if (!clingo_backend_rule(backend, true, head, 1, NULL, 0)) // choice rule, unfortunately not possible over strongly negated atoms are default negation
                EXIT_UNWANTED_STATE
        }

        // add clause
        clingo_literal_t body[clause.size()];
        // map literals
        for (int i = 0; i < (int)clause.size(); i++)
            body[i] = clause[i] > 0 ? -(clingo_literal_t)atoms[clause[i]] : atoms[-clause[i]]; // not that signs are reversed, so if all true then we have a problem.

        if (!clingo_backend_rule(backend, false, NULL, 0, body, clause.size())) //  choice rule unfortunately not possible over strongly negated atoms are default negation
            EXIT_UNWANTED_STATE
    }
}

bool ClingoSolver::propagate_clingo(clingo_propagate_control_t *control, const clingo_literal_t *lits, size_t n, ClingoSolver *s)
{
    // update current assignment based on lits
    for (int i = 0; i < (int)n; i++)
    {
        auto litsolver = lits[i];
        int lit = s->solverVariables2variables[abs(litsolver)];
        s->currentAssignment[lit] = litsolver > 0 ? truth_value_true : truth_value_false; // undo
    }

    s->propagate_control = control;
    if (!s->addStoredClauses())
        return true;
    s->propagate();
    return true;
}

bool ClingoSolver::undo_clingo(clingo_propagate_control_t *control, const clingo_literal_t *lits, size_t n, ClingoSolver *s)
{
    // remove assignemts from currentAssignment
    for (int i = 0; i < (int)n; i++)
    {
        int lit = s->solverVariables2variables[abs(lits[i])];
        s->currentAssignment[abs(lit)] = truth_value_unknown; // undo
    }
    return true;
}

bool ClingoSolver::check_clingo(clingo_propagate_control_t *control, ClingoSolver *s)
{
    const clingo_assignment_t *a = clingo_propagate_control_assignment(control);
    if (!clingo_assignment_is_total(a)) // clingo is buggy, sometimes not fully defined. Just to be sure
        return true;
    s->propagate_control = control;
    if (!s->addStoredClauses())
        return true;

    if (s->config.checkSolutionInProp)
    {
        s->check();
    }
    return true;
}

bool ClingoSolver::init_clingo(clingo_propagate_init_t *ctl, ClingoSolver *s)
{
    s->init(ctl);
    return true;
}
#include <typeinfo>
void ClingoSolver::init(clingo_propagate_init_t *ctl)
{
    solverVariables2variables = vector<clingo_literal_t>(variables2solverVariables.size(), 0);

    // TODO i think enough if I map all the observed literals to be able to update the current assignment
    for (auto l : config.observedVars)
    {
        clingo_atom_t edgeAtom = atoms[l];
        clingo_literal_t mappedLiteral;
        if (!clingo_propagate_init_solver_literal(ctl, edgeAtom, &mappedLiteral))
            EXIT_UNWANTED_STATE

        // printf("watch %d; %d %d; unmapped : %d\n", mappedEdges[i][j], i, j, edgeAtom);
        variables2solverVariables[l] = mappedLiteral;
        if (mappedLiteral < 0)
            EXIT_UNWANTED_STATE
        solverVariables2variables[mappedLiteral] = l;

        if (!clingo_propagate_init_add_watch(ctl, mappedLiteral))
            EXIT_UNWANTED_STATE
        if (!clingo_propagate_init_add_watch(ctl, -mappedLiteral))
            EXIT_UNWANTED_STATE
    }

    printf("Finished init\n");
    fflush(stdout);
}
