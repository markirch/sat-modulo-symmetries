#include "efx.hpp"

// #define EFX_VERSION1

#ifdef EFX_VERSION1
EFXPropagator::EFXPropagator(int frequency, int n, int p) : PartiallyDefinedGraphChecker(frequency), n(n), p(p)
{
    name = "EFX Propagator";

    // create encoding
    int nextVar = 1;
    int steps = p;
    edgeVars = vector<vector<int>>(n, vector<int>(n, 0));            // describes the directed graph
    vector<vector<int>> selectedPartition(steps, vector<int>(p, 0)); // indicates which partition is selected  in the s-th step
    selectedVertex = vector<vector<int>>(steps, vector<int>(n, 0));  // indicates which vertex is selected in the s-th step
    lenghtOfCycle = vector<int>(steps, 0);                           // indicates the length of the cycle
    vector<int> relevantStep(steps, 0);                              // indicates the relevant step (all steps < the length of the cycle)

    // create variables
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
                edgeVars[i][j] = nextVar++;

    for (int s = 0; s < steps; s++)
    {
        for (int i = 0; i < p; i++)
            selectedPartition[s][i] = nextVar++;
        for (int i = 0; i < n; i++)
            selectedVertex[s][i] = nextVar++;
        lenghtOfCycle[s] = nextVar++;
        relevantStep[s] = nextVar++;
    }

    cnf_t cnf;
    // in each step at least one partition is selected
    for (int s = 0; s < steps; s++)
        cnf.push_back(selectedPartition[s]);

    // each partition is selected at most once
    for (int i = 0; i < p; i++)
        for (int s = 0; s < steps; s++)
            for (int t = s + 1; t < steps; t++)
                cnf.push_back({-selectedPartition[s][i], -selectedPartition[t][i]});

    // in each step at least one vertex from the partition is selected
    int partitionSize = n / p;

    for (int s = 0; s < steps; s++)
    {
        for (int i = 0; i < p; i++)
        {
            clause_t clause = {-selectedPartition[s][i]};
            for (int j = partitionSize * i; j < partitionSize * (i + 1); j++)
            {
                clause.push_back(selectedVertex[s][j]);
            }
            cnf.push_back(clause);
        }
    }

    // if partition is not selected, no vertex from this partition is selected
    for (int s = 0; s < steps; s++)
        for (int i = 0; i < p; i++)
            for (int j = partitionSize * i; j < partitionSize * (i + 1); j++)
                cnf.push_back({selectedPartition[s][i], -selectedVertex[s][j]});

    // There must be an edge between the selected vertices in the cycle
    for (int s = 0; s < steps; s++)
        for (int s2 = s + 1; s2 < steps; s2++)
            cnf.push_back({-lenghtOfCycle[s2], relevantStep[s]}); // if cycle is larger then this must be considered

    for (int s = 0; s < steps - 1; s++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j)
                    cnf.push_back({-selectedVertex[s][i], -selectedVertex[s + 1][j], -relevantStep[s], edgeVars[i][j]});

    // The cycle must have a defined length
    cnf.push_back(lenghtOfCycle);
    cnf.push_back({-lenghtOfCycle[0]}); // cycles of length two should also be okay, because they are directed

    // The cycle must be closed
    for (int s = 1; s < steps; s++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j)
                    cnf.push_back({-lenghtOfCycle[s], -selectedVertex[s][i], -selectedVertex[0][j], edgeVars[i][j]}); // edge back to first selected vertex

    // load the clauses into the solver
    for (auto clause : cnf)
    {
        for (auto lit : clause)
            solver.add(lit);
        solver.add(0);
    }
}

void EFXPropagator::checkProperty(const adjacency_matrix_t &matrix)
{
    // Check if the graph is acyclic
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
            {
                if (matrix[i][j] == truth_value_true)
                    solver.assume(edgeVars[i][j]);
                else
                    solver.assume(-edgeVars[i][j]); // also undefined edges are assumed to be false because this doesn't introduce new cycles
            }

    if (solver.solve() != 10)
        return;

    // found cycle
    int length = 0;
    for (int i = 1; i < p; i++)
        if (solver.val(lenghtOfCycle[i]) > 0)
        {
            length = i + 1;
            break;
        }
    if (!length)
        EXIT_UNWANTED_STATE

    vector<int> cycle(length, -1);
    for (int s = 0; s < length; s++)
        for (int i = 0; i < n; i++)
            if (solver.val(selectedVertex[s][i]) > 0)
            {
                cycle[s] = i;
                break;
            }

    for (int i = 0; i < length; i++)
        if (cycle[i] == -1)
            EXIT_UNWANTED_STATE; // at least one vertex must be selected

    // print cycle
    // printf("Found cycle: ");
    // for (int i = 0; i < length; i++)
    //     printf("%d ", cycle[i]);
    // printf("\n");

    // forbid this directed cycle
    forbidden_graph_t forbiddenGraph;
    for (int i = 0; i < length - 1; i++)
        forbiddenGraph.push_back({truth_value_true, {cycle[i], cycle[i + 1]}});
    forbiddenGraph.push_back({truth_value_true, {cycle[length - 1], cycle[0]}});
    throw forbiddenGraph; // forbid this specific cycle
}

#else
// second version: only guess vertices in cycle and check if the incomming degree is 1 for all
EFXPropagator::EFXPropagator(int frequency, int n, int p) : PartiallyDefinedGraphChecker(frequency), n(n), p(p)
{
    name = "EFX Propagator";

    // create encoding
    int nextVar = 1;
    edgeVars = vector<vector<int>>(n, vector<int>(n, 0));       // describes the directed graph
    selectedVertices = vector<int>(n);                          // indicates which vertices are selected
    remainingEdges = vector<vector<int>>(n, vector<int>(n, 0)); // edges only between selected vertices

    // create variables
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
            {
                edgeVars[i][j] = nextVar++;
                remainingEdges[i][j] = nextVar++;
            }

    for (int i = 0; i < n; i++)
        selectedVertices[i] = nextVar++;

    cnf_t cnf;
    // in each step at least one vertex selected
    cnf.push_back(selectedVertices);

    // only one vertex of each partition is selected
    int partitionSize = n / p;
    for (int i = 0; i < p; i++)
    {
        for (int j = partitionSize * i; j < partitionSize * (i + 1); j++)
            for (int j2 = j + 1; j2 < partitionSize * (i + 1); j2++)
                cnf.push_back({-selectedVertices[j], -selectedVertices[j2]});
    }

    // set the remaining edges
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
            {
                cnf.push_back({-selectedVertices[i], -selectedVertices[j], -edgeVars[i][j], remainingEdges[i][j]}); // if selected by both and edge then present
                cnf.push_back({selectedVertices[i], -remainingEdges[i][j]});
                cnf.push_back({selectedVertices[j], -remainingEdges[i][j]});
                cnf.push_back({edgeVars[i][j], -remainingEdges[i][j]});
            }

    // ensure that the remaining degree is exactly 1 for selected vertices
    for (int i = 0; i < n; i++)
    {
        // if selected than at least one edge
        clause_t clause = {-selectedVertices[i]};
        for (int j = 0; j < n; j++)
            if (i != j)
                clause.push_back(remainingEdges[i][j]);
        cnf.push_back(clause); 
        // if selected than at least one incoming edge
        clause = {-selectedVertices[i]};
        for (int j = 0; j < n; j++)
            if (i != j)
                clause.push_back(remainingEdges[j][i]);
        cnf.push_back(clause);

        for (int j = 0; j < n; j++)
            if (i != j)
                for (int j2 = j + 1; j2 < n; j2++)
                    if (i != j2)
                        cnf.push_back({-remainingEdges[i][j], -remainingEdges[i][j2]}); // at most 1 independent whether selected or not
    }

    // load the clauses into the solver
    for (auto clause : cnf)
    {
        for (auto lit : clause)
            solver.add(lit);
        solver.add(0);
    }
}

void EFXPropagator::checkProperty(const adjacency_matrix_t &matrix)
{
    // Check if the graph is acyclic
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
            {
                if (matrix[i][j] == truth_value_true)
                    solver.assume(edgeVars[i][j]);
                else
                    solver.assume(-edgeVars[i][j]); // also undefined edges are assumed to be false because this doesn't introduce new cycles
            }

    if (solver.solve() != 10)
        return;

    vector<int> cycle;

    int startingVertex = -1;
    for (int i = 0; i < n; i++)
        if (solver.val(selectedVertices[i]) > 0)
        {
            startingVertex = i;
            break;
        }
    int nextVertex = startingVertex;
    do
    {
        // printf("next vertex: %d\n", nextVertex);
        cycle.push_back(nextVertex);
        for (int i = 0; i < n; i++)
            if (i != nextVertex && solver.val(remainingEdges[nextVertex][i]) > 0)
            {
                nextVertex = i;
                break;
            }
        if ((int) cycle.size() > n)
            EXIT_UNWANTED_STATE; // cycle is too long
    } while (nextVertex != startingVertex);

    // print cycle
    // printf("Found cycle version 2: ");
    // for (int i = 0; i < (int)cycle.size(); i++)
    //     printf("%d ", cycle[i]);
    // printf("\n");

    // forbid this directed cycle
    forbidden_graph_t forbiddenGraph;
    for (int i = 0; i < (int)cycle.size() - 1; i++)
        forbiddenGraph.push_back({truth_value_true, {cycle[i], cycle[i + 1]}});
    forbiddenGraph.push_back({truth_value_true, {cycle[cycle.size() - 1], cycle[0]}});
    throw forbiddenGraph; // forbid this specific cycle
}

#endif