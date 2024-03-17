/**
 * Used to construct decision trees with binary features
 * TODO generate by ChatGPT so check very carefully
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "decisionTree.hpp"

using namespace std;

// Function to calculate entropy
double entropy(const vector<sample_t> &samples)
{
    int num_samples = samples.size();
    if (num_samples == 0)
        return 0.0;

    int positive_count = 0;
    for (const auto &sample : samples)
    {
        if (sample.res)
        {
            positive_count++;
        }
    }

    double p_positive = static_cast<double>(positive_count) / num_samples;
    double p_negative = 1.0 - p_positive;

    if (p_positive == 0 || p_negative == 0)
    {
        return 0.0; // Avoid log(0)
    }

    return -p_positive * log2(p_positive) - p_negative * log2(p_negative);
}

// Function to find the best attribute to split on
int findBestAttribute(const vector<sample_t> &samples)
{
    int num_features = samples[0].features.size();
    double min_entropy = numeric_limits<double>::max();
    int best_attribute = -1;

    for (int i = 0; i < num_features; i++)
    {
        vector<sample_t> true_samples, false_samples;

        for (const auto &sample : samples)
        {
            if (sample.features[i])
            {
                true_samples.push_back(sample);
            }
            else
            {
                false_samples.push_back(sample);
            }
        }

        double weighted_entropy = (true_samples.size() * entropy(true_samples) +
                                   false_samples.size() * entropy(false_samples)) /
                                  samples.size();

        if (weighted_entropy < min_entropy)
        {
            min_entropy = weighted_entropy;
            best_attribute = i;
        }
    }

    return best_attribute;
}

// Recursive function to build the decision tree
Node *buildDecisionTree(const vector<sample_t> &samples)
{
    Node *node = new Node;

    // Check if all samples have the same result
    bool all_same_result = all_of(samples.begin(), samples.end(),
                                  [&](const sample_t &s)
                                  { return s.res == samples[0].res; });

    if (all_same_result)
    {
        node->is_leaf = true;
        node->label = samples[0].res;
    }
    else
    {
        node->is_leaf = false;
        node->split_attribute = findBestAttribute(samples);

        vector<sample_t> true_samples, false_samples;

        for (const auto &sample : samples)
        {
            if (sample.features[node->split_attribute])
            {
                true_samples.push_back(sample);
            }
            else
            {
                false_samples.push_back(sample);
            }
        }

        node->children.push_back(buildDecisionTree(true_samples));
        node->children.push_back(buildDecisionTree(false_samples));
    }

    return node;
}

void deleteDecisionTree(Node *node)
{
    if (node == nullptr)
    {
        return;
    }

    if (!node->is_leaf)
    {
        for (Node *child : node->children)
        {
            deleteDecisionTree(child);
        }
    }

    delete node;
}

void printDecisionTree(Node *node, int depth)
{
    if (node == nullptr)
        return;

    string indentation(depth, '\t');

    if (node->is_leaf)
    {
        cout << indentation << "Class: " << node->label << endl;
    }
    else
    {
        cout << indentation << "Attribute " << node->split_attribute << " at depth " << depth << endl;
        cout << indentation << "True branch:" << endl;
        printDecisionTree(node->children[0], depth + 1);
        cout << indentation << "False branch:" << endl;
        printDecisionTree(node->children[1], depth + 1);
    }
}

// Function to find all paths from the root to a leaf, including cases
void findAllPaths(Node *node, vector<pair<int, bool>> &path, vector<pair<vector<pair<int, bool>>, bool>> &paths)
{

    if (node->is_leaf)
    {
        paths.push_back({path, node->label}); // path + truth assignment at the end
    }
    else
    {
        path.push_back({node->split_attribute, true}); // TODO check if true or false
        findAllPaths(node->children[0], path, paths);
        path.pop_back();
        path.push_back({node->split_attribute, false});
        findAllPaths(node->children[1], path, paths);
        path.pop_back();
    }
}

/*

For testing

int main()
{
    // Sample data
    vector<sample_t> samples = {
        {{true, true, false, false}, true},
        {{true, true, true, false}, false},
        {{false, true, false, false}, true},
        {{false, false, true, false}, false},
        {{true, true, true, true}, true},
        {{false, false, false, true}, false}};

    // Build the decision tree
    Node *root = buildDecisionTree(samples);

    printDecisionTree(root);

    // Clean up memory (not shown in this simplified example)
    // You should implement a proper destructor for the Node structure.

    deleteDecisionTree(root);

    return 0;
}

*/
