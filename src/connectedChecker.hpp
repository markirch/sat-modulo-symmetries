#ifndef _FILTER_H_
#define _FILTER_H_

#include "useful.h"
#include "graphChecker.hpp"

class ConnectedChecker : public FullyDefinedGraphChecker
{
  void checkProperty(const adjacency_matrix_t &matrix);
public:
  ConnectedChecker()
  {
    name = "ConnectedChecker";
  }
};

#endif
