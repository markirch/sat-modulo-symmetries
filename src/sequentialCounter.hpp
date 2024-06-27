#ifndef SEQUENTIALCOUNTER_HPP
#define SEQUENTIALCOUNTER_HPP

#include "cadical.hpp"
#include "useful.h"

vector<int> sequentialCounter(CaDiCaL::Solver *solver, const vector<int> &vars, const int countUpTo, int &varCount, const int atLeast = -1, const int atMost = -1);


#endif
