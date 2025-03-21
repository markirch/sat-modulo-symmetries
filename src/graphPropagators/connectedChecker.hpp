#ifndef _CONNECTED_CHECKER_H_
#define _CONNECTED_CHECKER_H_

#include "../useful.h"
#include "../graphChecker.hpp"

class ConnectedChecker : public PartiallyDefinedGraphChecker
{
  void checkProperty(const adjacency_matrix_t &matrix);

public:
  ConnectedChecker()
  {
    name = "ConnectedChecker";
    this->onlyCheckFinal = true;
  }
};

class KConnectedChecker : public PartiallyDefinedGraphChecker
{
  int k;
  // ensure that graph is at least k connected
  void checkProperty(const adjacency_matrix_t &matrix);

public:
  KConnectedChecker(int k) : k(k)
  {
    name = "KConnectedChecker";
    this->onlyCheckFinal = true;
  }
};

#endif
