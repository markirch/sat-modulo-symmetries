#include "autcountChecker.hpp"
#include "../graph.hpp"

using std::make_pair;

void AutcountChecker::checkProperty(const adjacency_matrix_t &buggy_matrix, const vector<int> &model, int &nextFreeVariable)
{
    (void)buggy_matrix;     // Ignore the buggy matrix from SMS framework
    (void)nextFreeVariable; // Unused parameter
    
    // Workaround for SMS bug: Reconstruct correct matrix from model instead of using buggy_matrix
    // The SMS framework passes an incomplete matrix via currentAssignment, but the complete
    // model is available in the model vector. We reconstruct the correct matrix from it.
    
    int n = counter.getVertices(); // Get vertex count from counter
    GraphHandler temp_handler(n, false); // Assume undirected graphs
    
    // Convert model vector to truth_value_t assignment format
    // Note: We use model.size() but cap it to reasonable values to avoid SMS variable corruption
    size_t max_vars = std::min(model.size(), size_t(10000)); // Safety cap
    vector<truth_value_t> complete_assignment(max_vars, truth_value_unknown);
    
    for (size_t i = 1; i < max_vars; ++i) {
        if (model[i] > 0) {
            complete_assignment[abs(model[i])] = truth_value_true;
        } else {
            complete_assignment[abs(model[i])] = truth_value_false;
        }
    }
    
    // Reconstruct the correct, complete adjacency matrix
    adjacency_matrix_t correct_matrix = temp_handler.assignment2graph(complete_assignment);
    
    // Check if graph has at least minAutomorphisms using the correct matrix
    try {
        bool hasEnough = counter.hasAtLeastKAutomorphisms(correct_matrix, minAutomorphisms);
        
        if (!hasEnough) {
            // Graph doesn't have enough automorphisms - learn a clause
            // Create a blocking clause from the current graph assignment
            clause_t blocking_clause;
            
            // Get actual edge variables from the GraphHandler
            vector<vector<lit_t>> edge_vars = temp_handler.getEdgeVariables();
            for (int v1 = 0; v1 < n; ++v1) {
                for (int v2 = v1 + 1; v2 < n; ++v2) {
                    lit_t edge_var = edge_vars[v1][v2];
                    if (edge_var == 0) continue; // Skip if no variable assigned
                    
                    if (correct_matrix[v1][v2] == truth_value_true) {
                        blocking_clause.push_back(-edge_var);  // Negate true edges
                    } else if (correct_matrix[v1][v2] == truth_value_false) {
                        blocking_clause.push_back(edge_var);   // Negate false edges
                    }
                }
            }
            
            vector<clause_t> clauses;
            clauses.push_back(blocking_clause);
            throw clauses;
        }
    } catch (const std::runtime_error& e) {
        // VF2 computation failed - treat as not having enough automorphisms
        // Create the same blocking clause as above
        clause_t blocking_clause;
        vector<vector<lit_t>> edge_vars = temp_handler.getEdgeVariables();
        for (int v1 = 0; v1 < n; ++v1) {
            for (int v2 = v1 + 1; v2 < n; ++v2) {
                lit_t edge_var = edge_vars[v1][v2];
                if (edge_var == 0) continue;
                
                if (correct_matrix[v1][v2] == truth_value_true) {
                    blocking_clause.push_back(-edge_var);
                } else if (correct_matrix[v1][v2] == truth_value_false) {
                    blocking_clause.push_back(edge_var);
                }
            }
        }
        vector<clause_t> clauses;
        clauses.push_back(blocking_clause);
        throw clauses;
    }
}

void AutcountPDGChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    // For PDGs, check if ANY possible FDG completion can achieve minAutomorphisms
    // First check if the PDG has too many undefined edges for practical enumeration
    
    // Count undefined edges
    int undefined_count = 0;
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (matrix[i][j] == truth_value_unknown) {
                undefined_count++;
            }
        }
    }
    
    // PDG performance safeguard: if >15 undefined edges, be conservative and don't prune
    if (undefined_count > 15) {
        // Too many undefined edges for enumeration - assume PDG might achieve target
        // (Conservative approach: don't prune when we can't efficiently check)
        return;
    }
    
    try {
        bool canAchieveTarget = counter.hasAtLeastKAutomorphisms(matrix, minAutomorphisms);
        
        if (!canAchieveTarget) {
            // This PDG cannot possibly achieve the required automorphism count
            // Create forbidden_graph_t with all assigned edges as the conflict reason
            
            forbidden_graph_t forbidden_graph;
            int n = matrix.size();
            
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (matrix[i][j] != truth_value_unknown) {
                        // Add assigned edges (both true and false) to conflict reason
                        forbidden_graph.push_back(make_pair(matrix[i][j], make_pair(i, j)));
                    }
                }
            }
            
            throw forbidden_graph;
        }
    } catch (const std::runtime_error& e) {
        // Nauty computation failed on PDG (likely too many undefined edges)
        // Be conservative and create conflict to avoid missing violations
        
        forbidden_graph_t forbidden_graph;
        int n = matrix.size();
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (matrix[i][j] != truth_value_unknown) {
                    forbidden_graph.push_back(make_pair(matrix[i][j], make_pair(i, j)));
                }
            }
        }
        
        throw forbidden_graph;
    }
}