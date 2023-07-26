#ifndef PLANARITY_CHECKER_H
#define PLANARITY_CHECKER_H

#include "graphChecker.hpp"

// ensure that graph is planar
class PlanarityChecker : public PartiallyDefinedGraphChecker
{
public:
    PlanarityChecker(int frequency)
    {
        this->name = "PlanarityChecker";
        this->frequency = frequency;
    }
    void checkProperty(const adjacency_matrix_t &matrix);
};

// ensure that underlying graph of directed graph is planar
class DirectedPlanarityChecker : public PartiallyDefinedGraphChecker
{
public:
    DirectedPlanarityChecker(int frequency)
    {
        this->name = "DirectedPlanarityChecker";
        this->frequency = frequency;
    }
    void checkProperty(const adjacency_matrix_t &matrix);
};

/**ensure that thickness is two where the decomposition is represented by a directed graph where edges in both direction belong to the first graph and all single edges to the other one.
 By minimalty one can assume that unidirectional are from higher to lower vertices. (Must be encoded) */
class ThicknessTwoChecker : public PartiallyDefinedGraphChecker
{
public:
    ThicknessTwoChecker(int frequency)
    {
        this->name = "ThicknessTwoChecker";
        this->frequency = frequency;
    }
    void checkProperty(const adjacency_matrix_t &matrix);
};

class ThicknessTwoCheckerMulti : public PartiallyDefinedMultiGraphChecker
{
public:
    ThicknessTwoCheckerMulti(int frequency)
    {
        this->name = "ThicknessTwoCheckerMulti";
        this->frequency = frequency;
    }
    void checkProperty(const vector<adjacency_matrix_t> &matrices);
};



#endif