#include "autcount_maximizer.hpp"
#include "../sms.hpp"
#include "../graph.hpp"

using std::make_pair;

void AutcountMaximizerChecker::checkProperty(const adjacency_matrix_t &buggy_matrix, const vector<int> &model, int &nextFreeVariable)
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
    
    // Calculate exact automorphism count for this complete graph
    try {
        long long auts = counter.countExactAutomorphisms(correct_matrix);
        
        // Update solver's maximum and save graph if this is a new maximum or equal to current max
        solver.updateMaxAutomorphisms(static_cast<int>(auts), correct_matrix);
        
        // Always accept the graph (we're tracking maximum, not filtering)
        // The adaptive system works by updating the PDG pruning threshold
        
    } catch (const std::runtime_error& e) {
        // Nauty computation failed - still accept the graph but don't update maximum
        printf("⚠️  Automorphism computation failed for graph: %s\n", e.what());
    }
}

void AutcountMaximizerPDGChecker::checkProperty(const adjacency_matrix_t &matrix)
{
    // Get current maximum threshold from solver state
    int current_threshold = solver.getMaxAutomorphisms();
    
    // If threshold is still at initial value (0 or 1), don't prune anything yet
    if (current_threshold <= 1) {
        return; // Let search continue until we find some graphs to establish a baseline
    }
    
    // Count undefined edges for safety check
    int undefined_count = 0;
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (matrix[i][j] == truth_value_unknown) {
                undefined_count++;
            }
        }
    }
    
    // PDG performance safeguard: if too many undefined edges, be conservative and don't prune
    if (undefined_count > maxUndefinedEdges) {
        // Too many undefined edges for enumeration - assume PDG might achieve target
        // (Conservative approach: don't prune when we can't efficiently check)
        return;
    }
    
    try {
        // Check if this partial graph can possibly achieve the current maximum
        bool canAchieveThreshold = counter.hasAtLeastKAutomorphisms(matrix, current_threshold);
        
        if (!canAchieveThreshold) {
            // This PDG cannot possibly achieve the current maximum automorphism count
            // Create forbidden_graph_t with all assigned edges as the conflict reason
            
            forbidden_graph_t forbidden_graph;
            
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    if (matrix[i][j] != truth_value_unknown) {
                        // Add assigned edges (both true and false) to conflict reason
                        forbidden_graph.push_back(make_pair(matrix[i][j], make_pair(i, j)));
                    }
                }
            }
            
            // Throw conflict to prune this search branch
            throw forbidden_graph;
        }
    } catch (const std::runtime_error& e) {
        // Nauty computation failed on PDG (likely too many undefined edges)
        // Be conservative and don't prune (let search continue)
        return;
    }
}