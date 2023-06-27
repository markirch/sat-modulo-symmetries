#ifndef _CHECK_EMBEDDABILITY_H_
#define _CHECK_EMBEDDABILITY_H_

#include <vector>
#include "graphChecker.hpp"


bool isEmbeddable(std::vector<std::vector<bool>> G, int timeout = 2, int increment = 1);


// bool testEmebeddabilityKS(adjacency_matrix_t matrix)
// {
//   int n = matrix.size();
//   // printf("Test:\n");
//   // printAdjacencyMatrix(matrix);
//   vector<vector<bool>> adjacencyMatrix(n, vector<bool>(n, false));
//   for (int i = 0; i < n; i++)
//   {
//     for (int j = i + 1; j < n; j++)
//     {
//       if (matrix[i][j] == truth_value_true)
//         adjacencyMatrix[i][j] = adjacencyMatrix[j][i] = true;
//     }
//   }
//   clock_t start = clock();
//   auto x = isEmbeddable(adjacencyMatrix, 2, 1);
//   clock_t end = clock();
//   // printf("Time %f\n", ((double) end - start) / (CLOCKS_PER_SEC) );
//   // printf("ASFD %d\n", x ? 1 : 0);
//   return x;
// }

#endif