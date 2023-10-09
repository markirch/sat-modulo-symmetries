#ifndef _CONNECTED_CHECKER_H_
#define _CONNECTED_CHECKER_H_

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
