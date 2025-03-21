#include "embeddabilityCheck.hpp"

#include "z3++.h"
using namespace z3;

#include <queue>
#include <iostream>
#include <random>
#include <algorithm>
#include <string>


// Vector associated with a vertex in the graph
struct Vector
{
    std::vector<z3::expr> coordinates; // if empty then not assigned yet.
    int v1 = -1, v2 = -1;              //  potentially two adjacent vertices/orthogonal vectors given by index. -1 if not
};

// expr_vector
void crossproduct(Vector &v1, Vector &v2, std::vector<z3::expr> &coordinates)
{
    coordinates.push_back(v1.coordinates[1] * v2.coordinates[2] - v1.coordinates[2] * v2.coordinates[1]); // x
    coordinates.push_back(v1.coordinates[2] * v2.coordinates[0] - v1.coordinates[0] * v2.coordinates[2]); // y
    coordinates.push_back(v1.coordinates[0] * v2.coordinates[1] - v1.coordinates[1] * v2.coordinates[0]); // z
}

// check if we can assign to each vertex a vector. Timeout is given in milliseconds
check_result checkEmbeddibility(const std::vector<std::vector<bool>> &G, std::vector<Vector> &a, unsigned timeout, std::vector<int> &order, context &c)
{
    clock_t start = clock();
    solver s(c);
    clock_t end = clock();
    // printf("Time to create solver: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    int n = (int)G.size();

    printf("Size of graph to check: %d\n", n);

    for (int i = 0; i < n; i++)
    {
        // a[i].coordinates[0] = a[i].coordinates[0].simplify();
        // a[i].coordinates[1] = a[i].coordinates[1].simplify();
        // a[i].coordinates[2] = a[i].coordinates[2].simplify();
        // std::cout << "Crossproduct assignment: " << a[i].coordinates[0] << ";" << a[i].coordinates[1] << ";" << a[i].coordinates[2] << std::endl;
    }
    start = clock();

    for (auto v : order)
    {
        std::vector<z3::expr> cr;
        crossproduct(a[a[v].v1], a[a[v].v2], cr);
        s.add(a[v].coordinates[0] == cr[0]);
        s.add(a[v].coordinates[1] == cr[1]);
        s.add(a[v].coordinates[2] == cr[2]);
    }
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (G[i][j])
            {
                // printf("\torthogonal %d %d\n", i, j);
                // must be orthogonal
                s.add(a[i].coordinates[0] * a[j].coordinates[0] +
                          a[i].coordinates[1] * a[j].coordinates[1] +
                          a[i].coordinates[2] * a[j].coordinates[2] ==
                      0);
            }
        }
    end = clock();
    printf("Time for orth: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    start = clock();
    // vectors are not allowed to be collinear
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            // printf("\tcollinear %d %d\n", i, j);
            std::vector<z3::expr> cr;
            crossproduct(a[i], a[j], cr);
            s.add(cr[0] != 0 || cr[1] != 0 || cr[2] != 0); // if zero vector then the would be collinear
        }

    end = clock();
    printf("Time for collinear: %f\n", ((double)end - start) / CLOCKS_PER_SEC);

    // Create parameters object and set timeout to 10 seconds
    params p(c);
    p.set("timeout", timeout);
    // Set solver parameters
    s.set(p);
    end = clock();
    // printf("Time to set up: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    start = clock();
    auto res = s.check();
    end = clock();
    printf("Time to really solve: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    // std::cout << "Statistics: " << s.statistics() << std::endl;
    printf("Result for solving: %d\n", res);
    return res;
}

void assignCrossproductDependencies(const std::vector<std::vector<bool>> &G,
                                    std::vector<Vector> &assignment,
                                    std::deque<int> &expand,
                                    int n, std::vector<int> &order,
                                    context &ctx)
{

    while (expand.size() != 0)
    {
        int v = expand.front();
        expand.pop_front();

        for (int u = 0; u < n; ++u)
        {
            if (!assignment[u].coordinates.empty())
                continue; // already the coordinates are given

            if (!G[v][u])
                continue; // not adjacent so no impact

            if (assignment[u].v1 == -1)
                assignment[u].v1 = v; // one adjacent vector determined
            else if (assignment[u].v2 == -1)
            {
                assignment[u].v2 = v;
                order.push_back(u);
                auto x = "x" + std::to_string(u);
                auto y = "y" + std::to_string(u);
                auto z = "z" + std::to_string(u);
                assignment[u].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different
                expand.push_back(u);
                // order.push_back(u);
            }
            else
            {
                std::cout << "Should not happen" << std::endl;
            }
        }
    }
}

// remove low degree vertices and check connected components
std::vector<std::vector<bool>> preprocess(const std::vector<std::vector<bool>> &G)
{
    int n = (int)G.size();
    if (n < 10)
        return G;
    // printf("G %ld\n", G.size());
    // for (int i = 0; i < n; i++)
    //     printf("\t%ld\n", G[i].size());
    for (int i = 0; i < n; i++)
    {

        int degree;
        degree = 0;
        for (int j = 0; j < n; j++)
        {
            if (G[i][j])
                degree++;
        }

        if (degree <= 2)
        {
            std::vector<std::vector<bool>> Gnew(G.begin(), G.begin() + n - 1);
            for (int r = 0; r < n - 1; r++)
            {
                Gnew[r].pop_back();
            }

            if (i < n - 1)
            {

                for (int j = 0; j < n - 1; j++)
                {
                    if (j == i)
                        continue;
                    Gnew[i][j] = Gnew[i][j] = G[n - 1][j];
                }
            }
            return preprocess(Gnew);
        }
    }

    return G;
}

// produce a random assignment starting with a triangle and try to solve it with increasing timeout until we get a result
bool isEmbeddable(std::vector<std::vector<bool>> G, int timeout, int increment)
{
    G = preprocess(G);
    if (G.size() < 10)
        return true; // known to be embeddable
    printf("TEst\n");
    // auto r = test1(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }

    // r = test2(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }
    std::random_device rd;
    std::mt19937 g(rd());
    for (int round = 0; round < 100; round++) // at most 10 rounds
    {
        context ctx;
        int n = (int)G.size();
        // Randomized until solved
        // Find all triangles

        clock_t start = clock();
        std::vector<int> perm(n);
        for (int i = 0; i < n; i++)
            perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), g);

        std::vector<int> triangle;
        bool foundTriangle = false;
        // search for first triangle in the random permutation
        for (int i = 0; i < n && !foundTriangle; i++)
            for (int j = i + 1; j < n && !foundTriangle; j++)
            {

                if (!G[perm[i]][perm[j]])
                    continue; // definitely no triangle

                for (int k = j + 1; k < n; k++)
                {
                    if (!G[perm[i]][perm[k]] || !G[perm[j]][perm[k]])
                        continue; // no triangle
                    triangle = {perm[i], perm[j], perm[k]};
                    foundTriangle = true;
                }
            }
        // Gets one assignment
        std::vector<Vector> assignment(n, Vector());
        std::deque<int> expand;
        std::vector<int> order; // the order in which the cross products got assigned

        // Assign first triangle (if present)
        if (foundTriangle)
        {
            assignment[triangle[0]].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
            assignment[triangle[1]].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};
            assignment[triangle[2]].coordinates = {ctx.int_val(1), ctx.int_val(0), ctx.int_val(0)};
            expand = {triangle[0], triangle[1], triangle[2]};
        }
        else
        {
            // pick first edge in the random permutation
            bool foundEdge = false;
            for (int i = 0; i < n && !foundEdge; i++)
                for (int j = i + 1; j < n && !foundEdge; j++)
                {

                    if (!G[perm[i]][perm[j]])
                        continue;
                    foundEdge = true;
                    assignment[perm[i]].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
                    assignment[perm[j]].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};
                    expand = {perm[i], perm[j]};
                }
        }
        int numberOfFreeVectors = 0;
        assignCrossproductDependencies(G, assignment, expand, n, order, ctx);
        for (int i = 0; i < n; i++)
        {
            if (assignment[perm[i]].coordinates.empty())
            {
                numberOfFreeVectors++;
                auto x = "x" + std::to_string(perm[i]);
                auto y = "y" + std::to_string(perm[i]);
                auto z = "z" + std::to_string(perm[i]);

                assignment[perm[i]].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different
                expand = {perm[i]};
                assignCrossproductDependencies(G, assignment, expand, n, order, ctx);
            }
        }
        clock_t end = clock();
        printf("Number of free vectors: %d\n", numberOfFreeVectors);
        printf("Time for finding crossproduct cover: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
        start = clock();
        auto res = checkEmbeddibility(G, assignment, timeout, order, ctx);
        end = clock();
        printf("Real time for solving: %f\n\n", ((double)end - start) / CLOCKS_PER_SEC);

        if (res == sat)
        {
            return true;
        }
        else if (res == unsat)
        {
            return false;
        }
        // not solved; increase timeout and try other assignment
        timeout += increment;
    }

    return true;
}

/*


// Vector associated with a vertex in the graph
struct Vector
{
    std::vector<z3::expr> coordinates; // if empty then not assigned yet.
    int v1 = -1, v2 = -1;              //  potentially two adjacent vertices/orthogonal vectors given by index. -1 if not
};

// expr_vector
void crossproduct(Vector &v1, Vector &v2, std::vector<z3::expr> &coordinates)
{
    coordinates.push_back(v1.coordinates[1] * v2.coordinates[2] - v1.coordinates[2] * v2.coordinates[1]); // x
    coordinates.push_back(v1.coordinates[2] * v2.coordinates[0] - v1.coordinates[0] * v2.coordinates[2]); // y
    coordinates.push_back(v1.coordinates[0] * v2.coordinates[1] - v1.coordinates[1] * v2.coordinates[0]); // z
}

// check if we can assign to each vertex a vector. Timeout is given in milliseconds
check_result checkEmbeddibility(const std::vector<std::vector<bool>> &G, std::vector<Vector> &a, unsigned timeout, context &c)
{
    clock_t start = clock();
    solver s(c, solver::simple());
    clock_t end = clock();
    printf("Time to create solver: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    int n = (int)G.size();

    for (int i = 0; i < n; i++)
    {
        a[i].coordinates[0] = a[i].coordinates[0].simplify();
        a[i].coordinates[1] = a[i].coordinates[1].simplify();
        a[i].coordinates[2] = a[i].coordinates[2].simplify();
        // std::cout << "Crossproduct assignment: " << a[i].coordinates[0] << ";" << a[i].coordinates[1] << ";" << a[i].coordinates[2] << std::endl;
    }
    start = clock();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (G[i][j])
            {
                // printf("\torthogonal %d %d\n", i, j);
                // must be orthogonal
                s.add(a[i].coordinates[0] * a[j].coordinates[0] +
                          a[i].coordinates[1] * a[j].coordinates[1] +
                          a[i].coordinates[2] * a[j].coordinates[2] ==
                      0);
            }
        }
    end = clock();
    printf("Time for orth: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    start = clock();
    // vectors are not allowed to be collinear
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            // printf("\tcollinear %d %d\n", i, j);
            std::vector<z3::expr> cr;
            crossproduct(a[i], a[j], cr);
            s.add(cr[0] != 0 || cr[1] != 0 || cr[2] != 0); // if zero vector then the would be collinear
        }

    end = clock();
    printf("Time for collinear: %f\n", ((double)end - start) / CLOCKS_PER_SEC);

    // Create parameters object and set timeout to 10 seconds
    params p(c);
    p.set("timeout", timeout);
    // Set solver parameters
    s.set(p);
    end = clock();
    printf("Time to set up: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    start = clock();
    auto res = s.check();
    end = clock();
    printf("Time to really solve: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    std::cout << "Statistics: " << s.statistics() << std::endl;
    return res;
}

void assignCrossproductDependencies(const std::vector<std::vector<bool>> &G,
                                    std::vector<Vector> &assignment,
                                    std::deque<int> &expand,
                                    int n)
{

    while (expand.size() != 0)
    {
        int v = expand.front();
        expand.pop_front();

        for (int u = 0; u < n; ++u)
        {
            if (!assignment[u].coordinates.empty())
                continue; // already the coordinates are given

            if (!G[v][u])
                continue; // not adjacent so no impact

            if (assignment[u].v1 == -1)
                assignment[u].v1 = v; // one adjacent vector determined
            else if (assignment[u].v2 == -1)
            {
                assignment[u].v2 = v;
                crossproduct(assignment[assignment[u].v1], assignment[assignment[u].v2], assignment[u].coordinates);
                expand.push_back(u);
                // order.push_back(u);
            }
            else
            {
                std::cout << "Should not happen" << std::endl;
            }
        }
    }
}

// get assignment with as free as possible free vectors
check_result test1(const std::vector<std::vector<bool>> &G, int timeout)
{
    // First search for an edge then assing an additional vector; if only one free vector, then try to solve
    int n = (int)G.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            // printf("%d %d\n", i,j);
            if (!G[i][j])
                continue;

            for (int k = 0; k < n; k++)
            {
                if (k == i || k == j)
                    continue;

                // test
                context ctx;
                std::vector<Vector> assignment(n, Vector());
                assignment[i].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
                assignment[j].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};

                auto x = "x" + std::to_string(k);
                auto y = "y" + std::to_string(k);
                auto z = "z" + std::to_string(k);
                assignment[k].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different

                std::deque<int> expand = {i, j, k};
                assignCrossproductDependencies(G, assignment, expand, n);
                bool alreadyCovered = true;
                for (int l = 0; l < n; l++)
                {
                    if (assignment[l].coordinates.empty())
                    {
                        alreadyCovered = false;
                        break;
                    }
                }
                if (alreadyCovered)
                {
                    printf("Easy cover\n");
                    return checkEmbeddibility(G, assignment, timeout, ctx);
                }
            }
        }
    return unknown;
}

check_result test2(const std::vector<std::vector<bool>> &G, int timeout)
{
    // First search for an edge then assing an additional vector; if only one free vector, then try to solve
    int n = (int)G.size();
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            // printf("%d %d\n", i,j);
            if (!G[i][j])
                continue;

            for (int k = 0; k < n; k++)
            {
                if (k == i || k == j)
                    continue;

                for (int k2 = k + 1; k2 < n; k2++)
                {

                    // test
                    context ctx;
                    std::vector<Vector> assignment(n, Vector());
                    assignment[i].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
                    assignment[j].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};

                    auto x = "x" + std::to_string(k);
                    auto y = "y" + std::to_string(k);
                    auto z = "z" + std::to_string(k);
                    assignment[k].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different

                    x = "x" + std::to_string(k2);
                    y = "y" + std::to_string(k2);
                    z = "z" + std::to_string(k2);
                    assignment[k2].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different

                    std::deque<int> expand = {i, j, k, k2};
                    assignCrossproductDependencies(G, assignment, expand, n);
                    bool alreadyCovered = true;
                    for (int l = 0; l < n; l++)
                    {
                        if (assignment[l].coordinates.empty())
                        {
                            alreadyCovered = false;
                            break;
                        }
                    }
                    if (alreadyCovered)
                    {
                        printf("Easy cover2\n");
                        return checkEmbeddibility(G, assignment, timeout, ctx);
                    }
                }
            }
        }
    return unknown;
}

// remove low degree vertices and check connected components
std::vector<std::vector<bool>> preprocess(const std::vector<std::vector<bool>> &G)
{
    int n = (int)G.size();
    if (n < 10)
        return G;
    // printf("G %ld\n", G.size());
    // for (int i = 0; i < n; i++)
    //     printf("\t%ld\n", G[i].size());
    for (int i = 0; i < n; i++)
    {

        int degree;
        degree = 0;
        for (int j = 0; j < n; j++)
        {
            if (G[i][j])
                degree++;
        }

        if (degree <= 2)
        {
            std::vector<std::vector<bool>> Gnew(G.begin(), G.begin() + n - 1);
            for (int r = 0; r < n - 1; r++)
            {
                Gnew[r].pop_back();
            }

            if (i < n - 1)
            {

                for (int j = 0; j < n - 1; j++)
                {
                    if (j == i)
                        continue;
                    Gnew[i][j] = Gnew[i][j] = G[n - 1][j];
                }
            }
            return preprocess(Gnew);
        }
    }

    return G;
}

// produce a random assignment starting with a triangle and try to solve it with increasing timeout until we get a result
bool isEmbeddable(std::vector<std::vector<bool>> G, int timeout, int increment)
{
    G = preprocess(G);
    if (G.size() < 10)
        return true; // known to be embeddable
    // auto r = test1(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }

    // r = test2(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }
    std::random_device rd;
    std::mt19937 g(rd());
    for (int round = 0; round < 3; round++) // at most 10 rounds
    {
        context ctx;
        int n = (int)G.size();
        // Randomized until solved
        // Find all triangles

        clock_t start = clock();
        std::vector<int> perm(n);
        for (int i = 0; i < n; i++)
            perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), g);

        std::vector<int> triangle;
        bool foundTriangle = false;
        // search for first triangle in the random permutation
        for (int i = 0; i < n && !foundTriangle; i++)
            for (int j = i + 1; j < n && !foundTriangle; j++)
            {

                if (!G[perm[i]][perm[j]])
                    continue; // definitely no triangle

                for (int k = j + 1; k < n; k++)
                {
                    if (!G[perm[i]][perm[k]] || !G[perm[j]][perm[k]])
                        continue; // no triangle
                    triangle = {perm[i], perm[j], perm[k]};
                    foundTriangle = true;
                }
            }
        // Gets one assignment
        std::vector<Vector> assignment(n, Vector());
        std::deque<int> expand;

        // Assign first triangle (if present)
        if (foundTriangle)
        {
            assignment[triangle[0]].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
            assignment[triangle[1]].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};
            assignment[triangle[2]].coordinates = {ctx.int_val(1), ctx.int_val(0), ctx.int_val(0)};
            expand = {triangle[0], triangle[1], triangle[2]};
        }
        else
        {
            // pick first edge in the random permutation
            bool foundEdge = false;
            for (int i = 0; i < n && !foundEdge; i++)
                for (int j = i + 1; j < n && !foundEdge; j++)
                {

                    if (!G[perm[i]][perm[j]])
                        continue;
                    foundEdge = true;
                    assignment[perm[i]].coordinates = {ctx.int_val(0), ctx.int_val(0), ctx.int_val(1)};
                    assignment[perm[j]].coordinates = {ctx.int_val(0), ctx.int_val(1), ctx.int_val(0)};
                    expand = {perm[i], perm[j]};
                }
        }
        int numberOfFreeVectors = 0;
        assignCrossproductDependencies(G, assignment, expand, n);
        for (int i = 0; i < n; i++)
        {
            if (assignment[perm[i]].coordinates.empty())
            {
                numberOfFreeVectors++;
                auto x = "x" + std::to_string(perm[i]);
                auto y = "y" + std::to_string(perm[i]);
                auto z = "z" + std::to_string(perm[i]);

                assignment[perm[i]].coordinates = {ctx.real_const(x.c_str()), ctx.real_const(y.c_str()), ctx.real_const(z.c_str())}; // needs a name for variables but not necessary to be different
                expand = {perm[i]};
                assignCrossproductDependencies(G, assignment, expand, n);
            }
        }
        clock_t end = clock();
        printf("Number of free vectors: %d\n", numberOfFreeVectors);
        printf("Time for finding crossproduct cover: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
        start = clock();
        auto res = checkEmbeddibility(G, assignment, timeout, ctx);
        end = clock();
        printf("Real time for solving: %f\n", ((double)end - start) / CLOCKS_PER_SEC);

        if (res == sat)
        {
            return true;
        }
        else if (res == unsat)
        {
            return false;
        }
        // not solved; increase timeout and try other assignment
        timeout += increment;
    }

    return true;
}



#define MAX_VERTICES 24
bool alreadyIntialized = false;
context ctx;
solver s(ctx);

// Vector associated with a vertex in the graph
struct Vector
{
    int v1 = -1, v2 = -1; //  potentially two adjacent vertices/orthogonal vectors given by index. -1 if not
    bool assigned = false;
    int ort = 0; // if 1 then (0,0,1) if 2 then (0,1,0) if 3 then (1,0,0)
};

std::vector<expr_vector> coordinates; // coordinates which the crossproduct should be taken

void init()
{
    for (int i = 0; i < MAX_VERTICES; i++)
    {
        coordinates.push_back(expr_vector(ctx));

        auto x = "x" + std::to_string(i);
        auto y = "y" + std::to_string(i);
        auto z = "z" + std::to_string(i);
        coordinates[i].push_back(ctx.real_const(x.c_str()));
        coordinates[i].push_back(ctx.real_const(y.c_str()));
        coordinates[i].push_back(ctx.real_const(z.c_str()));
    }

    for (int i = 0; i < MAX_VERTICES; i++)
        for (int j = i + 1; j < MAX_VERTICES; j++)
        {
            // printf("\tcollinear %d %d\n", i, j);
            s.add(coordinates[i][1] * coordinates[j][2] - coordinates[i][2] * coordinates[j][1] != 0 ||
                  coordinates[i][2] * coordinates[j][0] - coordinates[i][0] * coordinates[j][2] != 0 ||
                  coordinates[i][0] * coordinates[j][1] - coordinates[i][1] * coordinates[j][0] != 0); // if zero vector then the would be collinear
        }
    alreadyIntialized = true;
}

// check if we can assign to each vertex a vector. Timeout is given in milliseconds
check_result checkEmbeddibility(const std::vector<std::vector<bool>> &G, std::vector<Vector> &a, unsigned timeout)
{
    if (!alreadyIntialized)
        init();

    clock_t start = clock();
    clock_t end = clock();
    // printf("Time to create solver: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    int n = (int)G.size();

    for (int i = 0; i < n; i++)
    {
        // a[i].coordinates[0] = a[i].coordinates[0].simplify();
        // a[i].coordinates[1] = a[i].coordinates[1].simplify();
        // a[i].coordinates[2] = a[i].coordinates[2].simplify();
        // std::cout << "Crossproduct assignment: " << a[i].coordinates[0] << ";" << a[i].coordinates[1] << ";" << a[i].coordinates[2] << std::endl;
    }
    start = clock();

    expr_vector assumptions(ctx);
    for (int i = 0; i < n; i++)
    {

        if (a[i].ort == 1)
        {
            assumptions.push_back(coordinates[i][0] == 0);
            assumptions.push_back(coordinates[i][1] == 0);
            assumptions.push_back(coordinates[i][2] == 1);
        }

        if (a[i].ort == 2)
        {
            assumptions.push_back(coordinates[i][0] == 0);
            assumptions.push_back(coordinates[i][1] == 1);
            assumptions.push_back(coordinates[i][2] == 0);
        }
        if (a[i].ort == 3)
        {
            assumptions.push_back(coordinates[i][0] == 1);
            assumptions.push_back(coordinates[i][1] == 0);
            assumptions.push_back(coordinates[i][2] == 0);
        }

        // croosproducts
        if (a[i].v1 != -1 && a[i].v2 != -1)
        {
            assumptions.push_back(coordinates[i][0] == coordinates[a[i].v1][1] * coordinates[a[i].v2][2] - coordinates[a[i].v1][2] * coordinates[a[i].v2][1]);
            assumptions.push_back(coordinates[i][1] == coordinates[a[i].v1][2] * coordinates[a[i].v2][0] - coordinates[a[i].v1][0] * coordinates[a[i].v2][2]);
            assumptions.push_back(coordinates[i][2] == coordinates[a[i].v1][0] * coordinates[a[i].v2][1] - coordinates[a[i].v1][1] * coordinates[a[i].v2][0]);
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
        {
            if (G[i][j])
            {
                // printf("\torthogonal %d %d\n", i, j);
                // must be orthogonal
                assumptions.push_back(coordinates[i][0] * coordinates[j][0] +
                                          coordinates[i][1] * coordinates[j][1] +
                                          coordinates[i][2] * coordinates[j][2] ==
                                      0);
            }
        }
    end = clock();
    printf("Time for orth + collinear: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    // Create parameters object and set timeout to 10 seconds
    params p(ctx);
    p.set("timeout", timeout);
    // Set solver parameters
    s.set(p);
    end = clock();
    // printf("Time to set up: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    start = clock();
    auto res = s.check(assumptions);
    end = clock();
    // printf("Time to really solve: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
    // std::cout << "Statistics: " << s.statistics() << std::endl;
    return res;
}

void assignCrossproductDependencies(const std::vector<std::vector<bool>> &G,
                                    std::vector<Vector> &assignment,
                                    std::deque<int> &expand,
                                    int n)
{

    while (expand.size() != 0)
    {
        int v = expand.front();
        expand.pop_front();

        for (int u = 0; u < n; ++u)
        {
            if (assignment[u].assigned)
                continue; // already assigned

            if (!G[v][u])
                continue; // not adjacent so no impact

            if (assignment[u].v1 == -1)
                assignment[u].v1 = v; // one adjacent vector determined
            else if (assignment[u].v2 == -1)
            {
                assignment[u].v2 = v;
                assignment[u].assigned = true;
                expand.push_back(u);
                // order.push_back(u);
            }
            else
            {
                std::cout << "Should not happen" << std::endl;
            }
        }
    }
}

// remove low degree vertices and check connected components
std::vector<std::vector<bool>> preprocess(const std::vector<std::vector<bool>> &G)
{
    int n = (int)G.size();
    if (n < 10)
        return G;
    // printf("G %ld\n", G.size());
    // for (int i = 0; i < n; i++)
    //     printf("\t%ld\n", G[i].size());
    for (int i = 0; i < n; i++)
    {

        int degree;
        degree = 0;
        for (int j = 0; j < n; j++)
        {
            if (G[i][j])
                degree++;
        }

        if (degree <= 2)
        {
            std::vector<std::vector<bool>> Gnew(G.begin(), G.begin() + n - 1);
            for (int r = 0; r < n - 1; r++)
            {
                Gnew[r].pop_back();
            }

            if (i < n - 1)
            {

                for (int j = 0; j < n - 1; j++)
                {
                    if (j == i)
                        continue;
                    Gnew[i][j] = Gnew[i][j] = G[n - 1][j];
                }
            }
            return preprocess(Gnew);
        }
    }

    return G;
}

// produce a random assignment starting with a triangle and try to solve it with increasing timeout until we get a result
bool isEmbeddable(std::vector<std::vector<bool>> G, int timeout, int increment)
{
    G = preprocess(G);
    if (G.size() < 10)
        return true; // known to be embeddable
    // auto r = test1(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }

    // r = test2(G, timeout);
    // if (r == sat)
    // {
    //     return true;
    // }
    // else if (r == unsat)
    // {
    //     return false;
    // }
    std::random_device rd;
    std::mt19937 g(rd());
    for (int round = 0; round < 3; round++) // at most 10 rounds
    {
        context ctx;
        int n = (int)G.size();
        // Randomized until solved
        // Find all triangles

        clock_t start = clock();
        std::vector<int> perm(n);
        for (int i = 0; i < n; i++)
            perm[i] = i;
        std::shuffle(perm.begin(), perm.end(), g);

        std::vector<int> triangle;
        bool foundTriangle = false;
        // search for first triangle in the random permutation
        for (int i = 0; i < n && !foundTriangle; i++)
            for (int j = i + 1; j < n && !foundTriangle; j++)
            {

                if (!G[perm[i]][perm[j]])
                    continue; // definitely no triangle

                for (int k = j + 1; k < n; k++)
                {
                    if (!G[perm[i]][perm[k]] || !G[perm[j]][perm[k]])
                        continue; // no triangle
                    triangle = {perm[i], perm[j], perm[k]};
                    foundTriangle = true;
                }
            }
        // Gets one assignment
        std::vector<Vector> assignment(n, Vector());
        std::deque<int> expand;
        std::vector<int> order; // the order in which the cross products got assigned

        // Assign first triangle (if present)
        if (foundTriangle)
        {
            assignment[triangle[0]].ort = 1;
            assignment[triangle[1]].ort = 2;
            assignment[triangle[2]].ort = 3;
            assignment[triangle[0]].assigned = true;
            assignment[triangle[1]].assigned = true;
            assignment[triangle[2]].assigned = true;
            expand = {triangle[0], triangle[1], triangle[2]};
        }
        else
        {
            // pick first edge in the random permutation
            bool foundEdge = false;
            for (int i = 0; i < n && !foundEdge; i++)
                for (int j = i + 1; j < n && !foundEdge; j++)
                {

                    if (!G[perm[i]][perm[j]])
                        continue;
                    foundEdge = true;
                    assignment[perm[i]].ort = 1;
                    assignment[perm[j]].ort = 2;
                    assignment[perm[i]].assigned = true;
                    assignment[perm[j]].assigned = true;
                    expand = {perm[i], perm[j]};
                }
        }
        int numberOfFreeVectors = 0;
        assignCrossproductDependencies(G, assignment, expand, n);
        for (int i = 0; i < n; i++)
        {
            if (!assignment[perm[i]].assigned)
            {
                numberOfFreeVectors++;
                assignment[perm[i]].assigned = true;
                expand = {perm[i]};
                assignCrossproductDependencies(G, assignment, expand, n);
            }
        }
        clock_t end = clock();
        printf("Number of free vectors: %d\n", numberOfFreeVectors);
        printf("Time for finding crossproduct cover: %f\n", ((double)end - start) / CLOCKS_PER_SEC);
        start = clock();
        auto res = checkEmbeddibility(G, assignment, timeout);
        end = clock();
        printf("Real time for solving: %f\n", ((double)end - start) / CLOCKS_PER_SEC);

        if (res == sat)
        {
            return true;
        }
        else if (res == unsat)
        {
            return false;
        }
        // not solved; increase timeout and try other assignment
        timeout += increment;
    }

    return true;
}

*/
