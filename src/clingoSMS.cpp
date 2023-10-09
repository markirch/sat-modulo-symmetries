#include "useful.h"
#include "clingoSMS.hpp"

// --------------------------------pure clingo part------------------------------------------

ClingoSolver::ClingoSolver(configSolver config, cnf_t &cnf) : GraphSolver(config)
{
    clingo_backend_t *backend;

    // create a control object and pass command line arguments
    // printf("Options %d\n", nOptionsClingo);
    // for (int i = 0; i < nOptionsClingo; i++)
    //     printf("%s\n", optionsClingo[i]);
    // if (!clingo_control_new(optionsClingo, nOptionsClingo, NULL, NULL, 20, &ctl) != 0)
    //     EXIT_UNWANTED_STATE

    // introduce as much variables as needed
    int highestVariable = 1;
    for (auto clause : cnf)
        for (auto lit : clause)
        {
            highestVariable = max(highestVariable, lit);
            highestVariable = max(highestVariable, -lit);
        }

    // get the backend
    if (!clingo_control_backend(ctl, &backend))
        EXIT_UNWANTED_STATE
    // prepare the backend for adding rules
    if (!clingo_backend_begin(backend))
        EXIT_UNWANTED_STATE
    // get atoms for each variable.

    atoms.push_back(0); // 0 is not an atom
    for (int i = 1; i <= highestVariable; i++)
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

void ClingoSolver::solve(vector<int> assumptions)
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

    if (!clingo_control_solve(ctl, clingo_solve_mode_yield, &mappedAssumptions[0], mappedAssumptions.size(), NULL, NULL, &handle))
        EXIT_UNWANTED_STATE
    if (!clingo_solve_handle_resume(handle))
        EXIT_UNWANTED_STATE
    if (!clingo_solve_handle_model(handle, &model))
        EXIT_UNWANTED_STATE
}

bool ClingoSolver::solve(vector<int> assumptions, int timeout)
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

    if (!clingo_control_solve(ctl, clingo_solve_mode_yield | clingo_solve_mode_async, &mappedAssumptions[0], mappedAssumptions.size(), NULL, NULL, &handle))
        EXIT_UNWANTED_STATE

    if (!clingo_solve_handle_resume(handle))
        EXIT_UNWANTED_STATE
    bool solvedSuccessfully;
    clingo_solve_handle_wait(handle, timeout, &solvedSuccessfully);
    return solvedSuccessfully;
}

adjacency_matrix_t ClingoSolver::getAdjacencyMatrix()
{
    adjacency_matrix_t matrix(vertices, std::vector<truth_value_t>(vertices, truth_value_false));
    const clingo_assignment_t *assignment = clingo_propagate_control_assignment(propagate_control);
#ifndef DIRECTED
    for (int i = 0; i < vertices; i++)
        for (int j = i + 1; j < vertices; j++)
        {
            clingo_truth_value_t value;
            if (!clingo_assignment_truth_value(assignment, mappedEdges[i][j], &value))
                EXIT_UNWANTED_STATE

            if (value == clingo_truth_value_free)
                matrix[j][i] = matrix[i][j] = truth_value_unknown;
            else if (value == clingo_truth_value_false)
                matrix[j][i] = matrix[i][j] = truth_value_false;
            else
                matrix[j][i] = matrix[i][j] = truth_value_true;
        }
    return matrix;
#else
    for (int i = 0; i < vertices; i++)
        for (int j = 0; j < vertices; j++)
        {
            if (i == j)
                continue;
            clingo_truth_value_t value;
            if (!clingo_assignment_truth_value(assignment, mappedEdges[i][j], &value))
                EXIT_UNWANTED_STATE

            if (value == clingo_truth_value_free)
                matrix[i][j] = truth_value_unknown;
            else if (value == clingo_truth_value_false)
                matrix[i][j] = truth_value_false;
            else
                matrix[i][j] = truth_value_true;
        }
    return matrix;
#endif
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

bool ClingoSolver::propagate_clingo(clingo_propagate_control_t *control, const clingo_literal_t *lits, size_t n, ClingoSolver *s)
{
    s->propagate_control = control;
    if (!s->addStoredClauses())
        return true;
    s->propagate();
    return true;
}

bool ClingoSolver::undo_clingo(clingo_propagate_control_t *control, const clingo_literal_t *, size_t, ClingoSolver *s)
{
    bool printBacktrack = false;
    if (printBacktrack)
    {
        EXIT_UNWANTED_STATE // control doesn't allow access while undo
            auto m = s->getAdjacencyMatrix();
        printPartiallyDefinedAdjacencyMatrix(m);
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
    s->check();
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
    if (!triangleVersion)
        if (!hyperedgeColoring)
            variables2solverVariables = vector<lit_t>(vertices * vertices + 1000, 0); // add 100 variables as threashold, not very rubost but should work
        else
            variables2solverVariables = vector<lit_t>(vertices * vertices + b_vertices[1] * b_vertices[1] + 1000, 0); // add 100 variables as threashold, not very rubost but should work
    else
        variables2solverVariables = vector<lit_t>(vertices * vertices * vertices + 1000, 0);

    mappedEdges = vector<vector<clingo_literal_t>>(vertices, vector<clingo_literal_t>(vertices));
    if (triangleVersion)
        mappedTriangles = vector<vector<vector<clingo_literal_t>>>(vertices, vector<vector<clingo_literal_t>>(vertices, vector<clingo_literal_t>(vertices)));
    if (hyperedgeColoring)
        mappedEdgesIntersectionGraph = vector<vector<clingo_literal_t>>(b_vertices[1], vector<clingo_literal_t>(b_vertices[1]));

#ifndef DIRECTED
    // map literals and add watches
    for (int i = 0; i < vertices; i++)
        for (int j = i + 1; j < vertices; j++)
        {
            clingo_atom_t edgeAtom = atoms[edges[i][j]];
            if (!clingo_propagate_init_solver_literal(ctl, edgeAtom, &(mappedEdges[i][j])))
                EXIT_UNWANTED_STATE

            mappedEdges[j][i] = mappedEdges[i][j];
            // printf("watch %d; %d %d; unmapped : %d\n", mappedEdges[i][j], i, j, edgeAtom);
            variables2solverVariables[edges[i][j]] = mappedEdges[i][j];

            if (!clingo_propagate_init_add_watch(ctl, mappedEdges[i][j]))
                EXIT_UNWANTED_STATE
            if (!clingo_propagate_init_add_watch(ctl, -mappedEdges[i][j]))
                EXIT_UNWANTED_STATE

            if (i < fixedSubgraphSize && j < fixedSubgraphSize)
                variablesOfSubgraphMapped.push_back(mappedEdges[i][j]);
        }

#else
    // map literals and add watches
    for (int i = 0; i < vertices; i++)
        for (int j = 0; j < vertices; j++)
        {
            if (i == j)
                continue;
            clingo_atom_t edgeAtom = atoms[edges[i][j]];
            if (!clingo_propagate_init_solver_literal(ctl, edgeAtom, &(mappedEdges[i][j])))
                EXIT_UNWANTED_STATE

            mappedEdges[i][j];
            // printf("watch %d; %d %d; unmapped : %d\n", mappedEdges[i][j], i, j, edgeAtom);
            variables2solverVariables[edges[i][j]] = mappedEdges[i][j];

            if (!clingo_propagate_init_add_watch(ctl, mappedEdges[i][j]))
                EXIT_UNWANTED_STATE
            if (!clingo_propagate_init_add_watch(ctl, -mappedEdges[i][j]))
                EXIT_UNWANTED_STATE

            if (i < fixedSubgraphSize && j < fixedSubgraphSize)
                variablesOfSubgraphMapped.push_back(mappedEdges[i][j]);
        }

#endif

    if (triangleVersion)
    {
        for (int i = 0; i < vertices; i++)
            for (int j = i + 1; j < vertices; j++)
                for (int k = j + 1; k < vertices; k++)
                {
                    clingo_atom_t triangleAtom = atoms[triangles[i][j][k]];
                    if (!clingo_propagate_init_solver_literal(ctl, triangleAtom, &(mappedTriangles[i][j][k])))
                        EXIT_UNWANTED_STATE

                    mappedTriangles[i][k][j] = mappedTriangles[j][i][k] = mappedTriangles[j][k][i] = mappedTriangles[k][i][j] = mappedTriangles[k][j][i] = mappedTriangles[i][j][k]; // assign for the rest
                    // printf("watch %d; %d %d; unmapped : %d\n", mappedEdges[i][j], i, j, triangleAtom);
                    variables2solverVariables[triangles[i][j][k]] = mappedTriangles[i][j][k];

                    if (!clingo_propagate_init_add_watch(ctl, mappedTriangles[i][j][k]))
                        EXIT_UNWANTED_STATE
                    if (!clingo_propagate_init_add_watch(ctl, -mappedTriangles[i][j][k]))
                        EXIT_UNWANTED_STATE

                    // if (i < fixedSubgraphSize && j < fixedSubgraphSize)
                    //     variablesOfSubgraphMapped.push_back(mappedEdges[i][j]);
                }
    }

    if (hyperedgeColoring)
    {
        // map literals and add watches
        for (int i = 0; i < b_vertices[1]; i++)
            for (int j = i + 1; j < b_vertices[1]; j++)
            {
                clingo_atom_t edgeAtom = atoms[edges_intersection_graph[i][j]];
                if (!clingo_propagate_init_solver_literal(ctl, edgeAtom, &(mappedEdgesIntersectionGraph[i][j])))
                    EXIT_UNWANTED_STATE

                mappedEdgesIntersectionGraph[j][i] = mappedEdgesIntersectionGraph[i][j];
                // printf("watch %d; %d %d; unmapped : %d\n", mappedEdges[i][j], i, j, edgeAtom);
                variables2solverVariables[edges_intersection_graph[i][j]] = mappedEdgesIntersectionGraph[i][j];

                if (!clingo_propagate_init_add_watch(ctl, mappedEdgesIntersectionGraph[i][j]))
                    EXIT_UNWANTED_STATE
                if (!clingo_propagate_init_add_watch(ctl, -mappedEdgesIntersectionGraph[i][j]))
                    EXIT_UNWANTED_STATE

                // if (i < fixedSubgraphSize && j < fixedSubgraphSize)
                //     variablesOfSubgraphMapped.push_back(mappedEdgesIntersectionGraph[i][j]);
            }
    }

    printf("Finished init\n");
    fflush(stdout);
}
