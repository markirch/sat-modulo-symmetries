#include "useful.h"
#include <numeric>
#include <string>
#include <string.h>
#include <sstream>

using std::string;

void printAdjacencyMatrix(const adjacency_matrix_t &matrix, bool printFullMatrix)
{
  if (printFullMatrix)
    for (auto row : matrix)
    {
      for (auto value : row)
        printf("%d ", value);
      printf("\n");
    }

  printf("[");
  bool first = true;

#ifndef DIRECTED
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = i + 1; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_true)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("]\n");
#else
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = 0; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_true)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("]\n");
#endif
}

void printPartiallyDefinedAdjacencyMatrix(const adjacency_matrix_t &matrix)
{

#ifndef DIRECTED
  printf("["); // other one
  printf("[");
  bool first = true;
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = i + 1; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_true)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("],");

  printf("[");
  first = true;
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = i + 1; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_false)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("],");

  printf("[");
  first = true;
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = i + 1; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_unknown)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("]");
  printf("]\n"); // other one
#else
  printf("[");
  bool first = true;
  for (int i = 0; i < (int)matrix.size(); i++)
    for (int j = 0; j < (int)matrix.size(); j++)
      if (matrix[i][j] == truth_value_true)
      {
        if (first)
          printf("(%d,%d)", i, j);
        else
          printf(",(%d,%d)", i, j);
        first = false;
      }
  printf("]\n");
#endif
}

void printHypergraph(const adjacency_matrix_t &matrix, int nVertices[2])
{
  char aux[3] = {',', '[', '{'};
  int och = 1;
  for (int e = nVertices[0]; e < nVertices[0] + nVertices[1]; e++)
  {
    printf("%c", aux[och]);
    och = 0;
    int ich = 2;
    for (int v = 0; v < nVertices[0]; v++)
    {
      if (matrix[v][e] == truth_value_true)
      {
        printf("%c%d", aux[ich], v);
        ich = 0;
      }
    }
    printf("}");
  }
  printf("]\n");
}

adjacency_matrix_t getIntersectionMatrix(const adjacency_matrix_t &incidence_matrix, int b_vertices[2])
{
  adjacency_matrix_t matrix(b_vertices[1], vector<truth_value_t>(b_vertices[1]));

  for (int e = 0; e < b_vertices[1]; e++)
  {
    for (int f = e + 1; f < b_vertices[1]; f++)
    {
      matrix[e][f] = matrix[f][e] = truth_value_false;
      for (int v = 0; v < b_vertices[0]; v++)
      {
        if (incidence_matrix[v][e + b_vertices[0]] && incidence_matrix[v][f + b_vertices[0]])
        {
          matrix[e][f] = matrix[f][e] = truth_value_true;
          break;
        }
      }
    }
  }
  return matrix;
}

// reads a file and stores the cnf and the maximum variable encountered
void file2cnf(std::ifstream &file, cnf_t &cnf, int &maxVar)
{
  maxVar = 0;
  clause_t clause; // current clause
  string line;
  while (getline(file, line))
  {
    if (strncmp(line.c_str(), "c\t", 2) == 0)
      continue;
    if (line[0] == 'p') // skips first line
      continue;
    std::istringstream iss(line);

    string space_delimiter = " ";

    string lit;
    while (std::getline(iss, lit, ' '))
    {
      int l = stoi(lit);
      maxVar = std::max(maxVar, abs(l));
      if (l == 0)
      {
        cnf.push_back(clause);
        clause = vector<int>(); // reset
      }
      else
        clause.push_back(l);
    }
  }
}
