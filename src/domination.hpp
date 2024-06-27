#ifndef DOMINATION_HPP
#define DOMINATION_HPP

#include "useful.h"
#include "graphChecker.hpp"

// special propagator for being quasi k-connected. Gives the number of interfaces and the size of the interfaces (all equal size). It is assumed that all interfaces are at the beginning.
class QuasiKConnectedPropagator : public FullyDefinedGraphChecker
{
    int interfaceSize;
    int numInterfaces;

    // constructor
public:
    QuasiKConnectedPropagator(int interfaceSize, int numInterfaces) : FullyDefinedGraphChecker()
    {
        this->name = "QuasiKConnectedPropagator";
        this->interfaceSize = interfaceSize;
        this->numInterfaces = numInterfaces;
    };

    // check property
    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif
