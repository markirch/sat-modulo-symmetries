#include "useful.h"

// Define the sample structure
struct sample_t
{
    vector<bool> features;
    bool res;
};

// Define a node for the decision tree
struct Node
{
    bool is_leaf;
    bool label;
    int split_attribute;
    vector<Node *> children;
};

// Recursive function to build the decision tree
Node *buildDecisionTree(const vector<sample_t> &samples);

void deleteDecisionTree(Node *node);

void printDecisionTree(Node *node, int depth = 0);

void findAllPaths(Node *node, vector<pair<int, bool>> &path, vector<pair<vector<pair<int, bool>>, bool>> &paths);
