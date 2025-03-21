#ifndef GRAPHCHECKER
#define GRAPHCHECKER

#include "useful.h"
#include <time.h>
#include <iostream>
#include <iomanip>

using std::string;

class GraphChecker
{
protected:
    string name = "Set name of Implementations of checher";
    long long calls = 0;           // number of times the property was checked
    clock_t time = 0;              // total time for checking the property
    long numberOfAddedClauses = 0; // total number of added clauses
    long long limit_reached = 0;   // count calls which did not provide a conclusive answer

public:
    void printStats()
    {
        std::cout << "Statistics for " << name << ":" << std::endl;
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "\tCalls: " << calls << std::endl;
        std::cout << "\tTime in seconds: " << ((double)time) / CLOCKS_PER_SEC << std::endl;
        std::cout << "\tAdded clauses: " << numberOfAddedClauses << std::endl;
        if (limit_reached > 0)
        {
            std::cout << "\tLimit reached: " << limit_reached << std::endl;
        }
    }
};

#define DEFAULT_FREQ 30

/**
 * A base class for checking partially defined graphs.
 * In case the property is violated, the edges/non-edges involved in the problem must be thrown as an exception
 */
class PartiallyDefinedGraphChecker : public GraphChecker
{
private:
    int counter = 0; // number of calls since last check

public:
    virtual ~PartiallyDefinedGraphChecker() = default;
    bool checkFinal = true;      // force check for fully defined graphs
    bool onlyCheckFinal = false; // If true, skip checks of partially defined graphs
    int frequency = DEFAULT_FREQ;

    PartiallyDefinedGraphChecker(int frequency = DEFAULT_FREQ) : frequency(frequency) {};

    /**
     * @brief The implementation fo the property to check
     *
     * @param matrix Current adjacency matrix
     * @throws forbidden_graph_t if the graph is forbidden
     */
    virtual void checkProperty(const adjacency_matrix_t &matrix) = 0;

    /**
     * @brief Check partially defined graph and throw partially defined graph which should be forbidden
     *
     * @param matrix The current partially defined graph given as adjacency matrix
     * @param final True if matrix is fully defined
     */
    void check(const adjacency_matrix_t &matrix, bool final)
    {
        counter++;
        if (final && !checkFinal)
            return;
        if (!final && counter % frequency != 0)
            return;

        if (!final && onlyCheckFinal)
            return;

        calls++;
        counter = 0;
        auto start = clock();
        try
        {
            checkProperty(matrix);
        }
        catch (const forbidden_graph_t forbiddenGraph)
        {
            time += clock() - start;
            numberOfAddedClauses++;
            throw forbiddenGraph;
        }
        time += clock() - start;
    }
};

/**
 * Allows more complicated checks at the end of the search.
 * Instead of a forbidden graph, the clauses are thrown. This can contain arbitrary variables and also introduce new once.
 */
class ComplexFullyDefinedGraphChecker : public GraphChecker
{
public:
    virtual ~ComplexFullyDefinedGraphChecker() = default;
    /**
     * @brief Implement the property which should be checked
     *
     * @param matrix The graph to check for
     * @param model The model given by cadical (i.e., true if > 0 and false if < 0)
     * @param nextFreeVariable The lowest available variable not used in the encoding or any other propagator
     * @throws vector<clause_t> if the property is violated
     */
    virtual void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable) = 0;

    /**
     * @brief Check graph and model and throw clauses which should be added. TODO think about giving the graph object
     *
     * @param matrix The current partially defined graph given as adjacency matrix
     * @param model The model
     * @param nextFreeVariable The lowest variable not used yet. If used, then must be incremented
     */
    void check(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable)
    {

        calls++;
        auto start = clock();
        try
        {
            checkProperty(matrix, model, nextFreeVariable);
        }
        catch (const vector<clause_t> clauses)
        {
            time += clock() - start;
            numberOfAddedClauses += clauses.size();
            throw clauses;
        }
        time += clock() - start;
    }
};

/**
 * Allows more complex propagators during search.
 * Note that all thrown literals must be observed !!!!!
 *
 */
class ComplexPartiallyDefinedGraphChecker : public GraphChecker
{
public:
    int counter; // number of calls since last check
    int frequency;
    bool checkFinal = true; // force check for fully defined graphs

    vector<int> observedVariables; // all potential variables which might be present in any thrown clause

    virtual void checkProperty(const adjacency_matrix_t &matrix, const vector<truth_value_t> &currentAssignemnt) = 0;

    /**
     * @brief Check partially defined graph and throw partially defined graph which should be forbidden
     *
     * @param matrix The current partially defined graph given as adjacency matrix
     * @param currentAssignemnt Assignemnt of the observed variables
     * @param final True if matrix is fully defined
     */
    void check(const adjacency_matrix_t &matrix, const vector<truth_value_t> &currentAssignemnt, bool final)
    {
        if (observedVariables.empty())
        {
            printf("Warning: no observed variables set in complex propagator.");
            EXIT_UNWANTED_STATE
        }

        counter++;
        if (final && !checkFinal)
            return;
        if (counter % frequency != 0 && !final)
            return;
        calls++;
        counter = 0;
        auto start = clock();
        try
        {
            checkProperty(matrix, currentAssignemnt);
        }
        catch (const vector<clause_t> clauses)
        {
            time += clock() - start;
            numberOfAddedClauses += clauses.size();
            throw clauses;
        }
        time += clock() - start;
    }
};

class PartiallyDefinedMultiGraphChecker : public GraphChecker
{
public:
    int counter = 0; // number of calls since last check
    int frequency = 30;
    bool checkFinal = true; // force check for fully defined graphs

    /**
     * @brief Check partially defined graph and throw partially defined graph which should be forbidden
     *
     * @param matrix The current partially defined graph given as adjacency matrix
     * @param final True if matrix is fully defined
     */
    void check(const vector<adjacency_matrix_t> &matrices, bool final)
    {
        counter++;
        if (final && !checkFinal)
            return;
        if (counter % frequency != 0 && !final)
            return;
        calls++;
        counter = 0;
        auto start = clock();
        try
        {
            checkProperty(matrices);
        }
        catch (const vector<forbidden_graph_t> forbiddenGraphs)
        {
            time += clock() - start;
            numberOfAddedClauses++;
            throw forbiddenGraphs;
        }
        time += clock() - start;
    }
    virtual void checkProperty(const vector<adjacency_matrix_t> &matrices) = 0;
};

#endif
