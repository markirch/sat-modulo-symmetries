#include "subgraphIsomorphism.hpp"
#include "cadical.hpp"
#include <math.h>

// currently only if target is also fully defined and non-induced
std::pair<cnf_t, vector<vector<int>>> createEncoding(const adjacency_matrix_t &pattern, const int nt)
{
    // target will be given as assumptions

    int np = (int)pattern.size();

    int varCount = 1;
    vector<vector<int>> edgeVars(nt, vector<int>(nt));
    for (int i = 0; i < nt; i++)
        for (int j = i + 1; j < nt; j++)
            edgeVars[j][i] = edgeVars[i][j] = varCount++;

    vector<vector<int>> mapping(np, vector<int>(nt));
    for (int i = 0; i < np; i++)
        for (int j = 0; j < nt; j++)
            mapping[i][j] = varCount++;

    cnf_t cnf;
    // for each vertex in pattern, there is at least one vertex in target
    for (int i = 0; i < np; i++)
        cnf.push_back(mapping[i]);
    // injectivity
    for (int i = 0; i < np; i++)
        for (int j = i + 1; j < np; j++)
            for (int k = 0; k < nt; k++)
                cnf.push_back({-mapping[i][k], -mapping[j][k]});
    // at most one vertex in target for each vertex in pattern
    for (int i = 0; i < np; i++)
        for (int j = 0; j < nt; j++)
            for (int k = j + 1; k < nt; k++)
                cnf.push_back({-mapping[i][j], -mapping[i][k]});

    // edges must map to edges  TODO non induced case not handeled yet
    for (int i = 0; i < np; i++)
        for (int j = i + 1; j < np; j++)
            for (int k = 0; k < nt; k++)
                for (int l = k + 1; l < nt; l++)
                {
                    if (pattern[i][j] == truth_value_true)
                    {
                        cnf.push_back({-mapping[i][k], -mapping[j][l], edgeVars[k][l]});
                        cnf.push_back({-mapping[i][l], -mapping[j][k], edgeVars[k][l]});
                    }
                }

    // cardinality constraints for neighborhoods
    if (1)
        for (int i = 0; i < np; i++)
        {
            vector<int> isNeighbor(nt);
            for (int j = 0; j < nt; j++)
                isNeighbor[j] = varCount++;

            vector<int> neighborsInPattern;
            for (int j = 0; j < np; j++)
                if (pattern[i][j] == truth_value_true)
                    neighborsInPattern.push_back(j);

            // if isNeighbor[j] is true than on neighbor of i must be mapped to j
            for (int j = 0; j < nt; j++)
            {
                vector<int> clause = {-isNeighbor[j]};
                for (int k : neighborsInPattern)
                    clause.push_back(mapping[k][j]);
                cnf.push_back(clause);
            }

            // sequential counter for ensuring that at least deg(i) variable of isNeighbor are true
            int deg = (int)neighborsInPattern.size();
            int countTo = deg;
            vector<vector<int>> counterVars(nt, vector<int>(countTo));
            for (int j = 0; j < nt; j++)
                for (int k = 0; k < countTo; k++)
                    counterVars[j][k] = varCount++;

            // base case
            counterVars[0][0] = isNeighbor[0];
            for (int k = 1; k < countTo; k++)
                cnf.push_back({-counterVars[0][k]});

            for (int j = 1; j < nt; j++)
            {
                for (int k = 0; k < countTo; k++)
                {
                    cnf.push_back({isNeighbor[j], counterVars[j - 1][k], -counterVars[j][k]}); // if it isn't a neighbor and the previous counter isn't true, then the current counter must be false
                    if (k > 0)
                        cnf.push_back({counterVars[j - 1][k - 1], -counterVars[j][k]}); // if the previous counter - 1 is false, then the current counter must be false
                }
            }

            cnf.push_back({counterVars[nt - 1][countTo - 1]}); // at least deg(i) neighbors
        }
    return {cnf, edgeVars};
}

// Returns true when the k-vertex subgraph (with adjacency matrix M) contains the gth minimal unembeddable graph
// M is determined by the current assignment to the first k*(k-1)/2 variables
// If an unembeddable subgraph is found in M, P is set to the mapping from the rows of M to the rows of the lex greatest representation of the unembeddable adjacency matrix (and p is the inverse of P)
bool has_gub_subgraph(adjacency_matrix_t pattern, adjacency_matrix_t target)
{
    int k = (int) target.size();

    int pl[12]; // pl[k] contains the current list of possibilities for kth vertex (encoded bitwise)
    int pn[13]; // pn[k] contains the initial list of possibilities for kth vertex (encoded bitwise)
    pl[0] = (1 << k) - 1;
    pn[0] = (1 << k) - 1;
    int i = 0;

    vector<int> p(pattern.size());
    vector<int> P(target.size());

    while (1)
    {
        // If no possibilities for ith vertex then backtrack
        if (pl[i] == 0)
        {
            // Backtrack to vertex that has at least two possibilities
            while ((pl[i] & (pl[i] - 1)) == 0)
            {
                i--;
                if (i == -1)
                {
                    // No permutations produce a matrix containing the gth submatrix
                    return false;
                }
            }
            // Remove p[i] as a possibility from the ith vertex
            pl[i] = pl[i] & ~(1 << p[i]);
        }

        p[i] = log2(pl[i] & -pl[i]);      // Get index of rightmost high bit
        pn[i + 1] = pn[i] & ~(1 << p[i]); // List of possibilities for (i+1)th vertex

        // Determine if the permuted matrix p(M) is contains the gth submatrix
        bool result_known = false;
        for (int j = 0; j < i; j++)
        {
            if (pattern[j][i] == truth_value_false) //  !gub[g][i*(i-1)/2+j])
                continue;
            if (target[p[i]][p[j]] != truth_value_true)
            {
                // Permutation sends a non-edge to a gth submatrix edge; stop considering
                result_known = true;
                break;
            }
        }

        if (!result_known && i == (int)pattern.size() - 1)
        {
            // The complete gth submatrix found in p(M)
            for (int j = 0; j <= i; j++)
            {
                P[p[j]] = j;
            }
            return true;
        }
        if (!result_known)
        {
            // Result is unknown; need to define p[i] for another i
            i++;
            pl[i] = pn[i];
        }
        else
        {
            // Remove p[i] as a possibility from the ith vertex
            pl[i] = pl[i] & ~(1 << p[i]);
        }
    }
}

clock_t totalCadicalTime = 0;
clock_t totalSimpleTime = 0;
void ForbiddenSubgraphChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    int n = (int)matrix.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (matrix[i][j] == truth_value_true)
                forbiddenSubgraphSolver.assume(edgeVars[i][j]);
            else
                forbiddenSubgraphSolver.assume(-edgeVars[i][j]);

    // time precisely solver
    // clock_t start = clock();
    /*int res =*/ forbiddenSubgraphSolver.solve();

    printf("ERROR: ForbiddenSubgraphChecker using Cadical is not fully implemented yet\n");
    EXIT_UNWANTED_STATE // still work in progress
    // printf("result: %d\n", res);
    // printf("ForbiddenSubgraphChecker: %f\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    // totalCadicalTime += clock() - start;

    // start = clock();
    // adjacency_matrix_t matrix_copy = matrix;
    // // set all undefined edges to false
    // for (int i = 0; i < n; i++)
    //     for (int j = i + 1; j < n; j++)
    //         if (matrix[i][j] == truth_value_unknown)
    //             matrix_copy[i][j] = matrix_copy[j][i] = truth_value_false;
    // bool res2 = has_gub_subgraph(forbiddenSubgraphs[0], matrix);
    // printf("result: %d\n", res);
    // printf("Simple: %f\n", (double)(clock() - start) / CLOCKS_PER_SEC);
    // totalSimpleTime += clock() - start;

    // printf("Ratio: %f\n", (double)totalSimpleTime / totalCadicalTime);
}

void ForbiddenSubgraphChecker::load_graphs(std::ifstream &graphs, bool /*induced*/)
{
    string line;
    while (getline(graphs, line))
    {
        if (strncmp(line.c_str(), "c\t", 2) == 0)
            continue;

        vector<std::pair<int, int>> edges;

        std::istringstream iss(line);
        clause_t clause;
        string space_delimiter = " ";

        string lit;
        std::getline(iss, lit, ' ');
        int n = stoi(lit); // first integer gives number of vertices the remaining the edges

        adjacency_matrix_t m(n, vector<truth_value_t>(n, truth_value_false));

        while (std::getline(iss, lit, ' '))
        {
            int v1 = stoi(lit);
            std::getline(iss, lit, ' ');
            int v2 = stoi(lit);
            m[v1][v2] = m[v2][v1] = truth_value_true;
        }
        forbiddenSubgraphs.push_back(m);
    }
}
