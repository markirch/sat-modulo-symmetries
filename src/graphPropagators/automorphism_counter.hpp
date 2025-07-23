#ifndef AUTOMORPHISM_COUNTER_HPP
#define AUTOMORPHISM_COUNTER_HPP

#include "../useful.h"
#include <igraph/igraph.h>
#include <vector>
#include <unordered_map>

class AutomorphismCounter {
private:
    int vertices;
    bool directed;
    
    // Early termination context for igraph callbacks
    struct EarlyTerminationContext {
        int target_count;
        int current_count;  
        bool should_stop;
        const adjacency_matrix_t* matrix1;
        const adjacency_matrix_t* matrix2;
    };
    
    // Pure C++ helper methods
    int countTrueNeighbors(const adjacency_matrix_t& matrix, vertex_t v) const;
    int countZeroNeighborVertices(const adjacency_matrix_t& matrix) const;
    long long computeFactorial(int n) const;
    
    // SMS adjacency matrix ↔ igraph conversion
    igraph_t convertToIgraph(const adjacency_matrix_t& matrix);
    void destroyIgraph(igraph_t& graph);
    
    // VF2 with early termination implementation
    bool vf2EarlyTermination(const adjacency_matrix_t& matrix, int k);
    
    
    // Static callbacks for igraph (C API requirement)
    static igraph_error_t early_termination_callback(
        const igraph_vector_int_t *map12,
        const igraph_vector_int_t *map21, 
        void *arg);
        
    static bool edge_compatibility_callback(
        const igraph_t *graph1,
        const igraph_t *graph2,
        igraph_integer_t vertex1,
        igraph_integer_t vertex2,
        void *arg);

public:
    explicit AutomorphismCounter(int v, bool dir = false);
    ~AutomorphismCounter();
    
    // Main interface methods
    bool hasAtLeastKAutomorphisms(const adjacency_matrix_t& matrix, int k);
    long long countExactAutomorphisms(const adjacency_matrix_t& matrix);
    
    // Getter methods
    int getVertices() const { return vertices; }
    
    // Pure C++ bypass (no external dependencies)
    bool simpleBypass(const adjacency_matrix_t& matrix, int k) const;
};

// Pure C++ SMS edge compatibility helper (outside class)
static bool areEdgesCompatible(truth_value_t edge1, truth_value_t edge2) {
    // SMS PDG compatibility rules:
    // Allowed: 1↔1, 0↔0, *↔*, 1↔*, 0↔*, *↔1, *↔0  
    // Forbidden: 1↔0, 0↔1
    
    if (edge1 == truth_value_true && edge2 == truth_value_false) return false;
    if (edge1 == truth_value_false && edge2 == truth_value_true) return false;
    return true;
}

#endif // AUTOMORPHISM_COUNTER_HPP