#ifndef UNSAT_CHECKER_H
#define UNSAT_CHECHER_H
#include "cadical.hpp"
#include "graphChecker.hpp"
#include "useful.h"

class UNSATChecker : public ComplexFullyDefinedGraphChecker
{   
    int nof_var;
    int nof_cls;
    vector<vector<int>> *edge_stats;

public:
    UNSATChecker(int n_var, int n_cls)
    {
        this->nof_var = n_var;
        this->nof_cls = n_cls;
        this->edge_stats = edge_stats;
        this->name = "UNSATChecker";
        addsOnlyObservedLiterals = false;
    }
    void checkProperty(const adjacency_matrix_t &matrix);
};

#endif