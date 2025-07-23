#include "automorphism_counter.hpp"
#include <iostream>
#include <stdexcept>
#include <limits>

AutomorphismCounter::AutomorphismCounter(int v, bool dir) 
    : vertices(v), directed(dir) {
    if (vertices <= 0) {
        throw std::invalid_argument("Vertex count must be positive");
    }
}

AutomorphismCounter::~AutomorphismCounter() {
    // Clean up resources if needed
}

// Simple zero 1-neighbor bypass implementation
bool AutomorphismCounter::simpleBypass(const adjacency_matrix_t& matrix, int k) const {
    // Trivial case
    if (k <= 1) return true;
    
    // Count vertices with zero confirmed edges (truth_value_true)
    int zero_neighbor_count = 0;
    
    for (vertex_t v = 0; v < vertices; v++) {
        int true_neighbors = countTrueNeighbors(matrix, v);
        if (true_neighbors == 0) {
            zero_neighbor_count++;
        }
    }
    
    // Check if factorial(count) >= k
    long long factorial_bound = computeFactorial(zero_neighbor_count);
    return factorial_bound >= k;
}

int AutomorphismCounter::countTrueNeighbors(const adjacency_matrix_t& matrix, vertex_t v) const {
    int count = 0;
    for (vertex_t u = 0; u < vertices; u++) {
        if (v != u && matrix[v][u] == truth_value_true) {
            count++;
        }
    }
    return count;
}

int AutomorphismCounter::countZeroNeighborVertices(const adjacency_matrix_t& matrix) const {
    int count = 0;
    for (vertex_t v = 0; v < vertices; v++) {
        if (countTrueNeighbors(matrix, v) == 0) {
            count++;
        }
    }
    return count;
}

long long AutomorphismCounter::computeFactorial(int n) const {
    const long long MAX_FACTORIAL = 1000000000LL; // Overflow protection
    
    if (n <= 0) return 1;  
    if (n > 20) return MAX_FACTORIAL; // 21! > 2^64
    
    long long result = 1;
    for (int i = 2; i <= n; i++) {
        if (result > MAX_FACTORIAL / i) return MAX_FACTORIAL;
        result *= i;
    }
    return result;
}

// SMS adjacency matrix to igraph conversion
igraph_t AutomorphismCounter::convertToIgraph(const adjacency_matrix_t& matrix) {
    igraph_t graph;
    igraph_vector_int_t edges;
    
    igraph_vector_int_init(&edges, 0);
    
    // Check if this is a PDG (has unknown edges)
    bool hasUnknownEdges = false;
    for (vertex_t i = 0; i < vertices && !hasUnknownEdges; i++) {
        for (vertex_t j = 0; j < vertices && !hasUnknownEdges; j++) {
            if (matrix[i][j] == truth_value_unknown) {
                hasUnknownEdges = true;
            }
        }
    }
    
    if (hasUnknownEdges) {
        // For PDGs: create complete clique with ALL possible edges
        // VF2 edge compatibility callback will handle constraint checking
        for (vertex_t i = 0; i < vertices; i++) {
            for (vertex_t j = (directed ? 0 : i + 1); j < vertices; j++) {
                if (i != j) {
                    igraph_vector_int_push_back(&edges, i);
                    igraph_vector_int_push_back(&edges, j);
                }
            }
        }
    } else {
        // For fully-defined graphs: create actual graph structure
        for (vertex_t i = 0; i < vertices; i++) {
            for (vertex_t j = (directed ? 0 : i + 1); j < vertices; j++) {
                if (i != j && matrix[i][j] == truth_value_true) {
                    igraph_vector_int_push_back(&edges, i);
                    igraph_vector_int_push_back(&edges, j);
                }
            }
        }
    }
    
    igraph_create(&graph, &edges, vertices, directed);
    igraph_vector_int_destroy(&edges);
    
    return graph;
}

void AutomorphismCounter::destroyIgraph(igraph_t& graph) {
    igraph_destroy(&graph);
}

// Static callbacks for igraph  
igraph_error_t AutomorphismCounter::early_termination_callback(
    const igraph_vector_int_t *map12,
    const igraph_vector_int_t *map21, 
    void *arg) {
    
    (void)map12; // Suppress unused parameter warning
    (void)map21; // Suppress unused parameter warning
    
    EarlyTerminationContext* ctx = static_cast<EarlyTerminationContext*>(arg);
    ctx->current_count++;
    
    if (ctx->current_count >= ctx->target_count) {
        ctx->should_stop = true;
        return IGRAPH_STOP; // Early termination
    }
    
    return IGRAPH_SUCCESS;
}

bool AutomorphismCounter::edge_compatibility_callback(
    const igraph_t *graph1, const igraph_t *graph2,
    igraph_integer_t vertex1, igraph_integer_t vertex2, void *arg) {
    
    (void)graph1; // We'll use the stored matrices instead
    (void)graph2;
    
    EarlyTerminationContext* ctx = static_cast<EarlyTerminationContext*>(arg);
    
    // Get truth values from SMS adjacency matrix using vertex IDs
    truth_value_t val1 = (*ctx->matrix1)[vertex1][vertex2];  
    truth_value_t val2 = (*ctx->matrix2)[vertex1][vertex2];  // Same for automorphism
    
    return areEdgesCompatible(val1, val2);
}

// VF2 early termination context for callbacks
struct VF2EarlyTerminationContext {
    long long target_count;
    long long current_count;
    const adjacency_matrix_t* matrix;
    bool should_stop;
};

// VF2 edge compatibility callback - implements SMS PDG edge compatibility
extern "C" igraph_bool_t vf2_edge_compatibility_callback(
    const igraph_t* graph1, const igraph_t* graph2,
    igraph_integer_t edge1, igraph_integer_t edge2, void* arg) {
    
    VF2EarlyTerminationContext* ctx = static_cast<VF2EarlyTerminationContext*>(arg);
    
    igraph_integer_t u1, v1, u2, v2;
    
    // Get vertices corresponding to edge1 in graph1
    igraph_edge(graph1, edge1, &u1, &v1);
    // Get vertices corresponding to edge2 in graph2
    igraph_edge(graph2, edge2, &u2, &v2);
    
    // Retrieve truth values from the original adjacency matrix using the vertex indices.
    // For automorphism, both edges come from the same matrix.
    truth_value_t val1 = (*ctx->matrix)[u1][v1];
    truth_value_t val2 = (*ctx->matrix)[u2][v2];
    
    return areEdgesCompatible(val1, val2);
}

// VF2 automorphism found callback - counts and implements early termination
extern "C" igraph_error_t vf2_automorphism_found_callback(
    const igraph_vector_int_t* map12, const igraph_vector_int_t* map21, void* arg) {
    
    (void)map12; (void)map21; // We just count for now
    
    VF2EarlyTerminationContext* ctx = static_cast<VF2EarlyTerminationContext*>(arg);
    ctx->current_count++;
    
    // Early termination: stop if we've found enough automorphisms
    if (ctx->current_count >= ctx->target_count) {
        ctx->should_stop = true;
        return IGRAPH_STOP; // Stop search
    }
    
    return IGRAPH_SUCCESS; // Continue search
}

// VF2 with early termination implementation  
bool AutomorphismCounter::vf2EarlyTermination(const adjacency_matrix_t& matrix, int k) {
    // Convert SMS matrix to igraph
    igraph_t graph = convertToIgraph(matrix);
    
    // Set up early termination context
    VF2EarlyTerminationContext ctx;
    ctx.target_count = k;
    ctx.current_count = 0;
    ctx.matrix = &matrix;
    ctx.should_stop = false;
    
    // Use VF2 for automorphism counting with early termination
    igraph_error_t result = igraph_isomorphic_function_vf2(
        &graph, &graph,  // Same graph (automorphism)
        nullptr, nullptr,  // No vertex colors
        nullptr, nullptr,  // No edge colors  
        nullptr, nullptr,  // No mapping output needed
        vf2_automorphism_found_callback,  // Count automorphisms
        nullptr,  // No vertex compatibility needed
        vf2_edge_compatibility_callback,  // SMS edge compatibility
        &ctx  // Our context
    );
    
    destroyIgraph(graph);
    
    // IGRAPH_STOP is a valid return code when early termination is triggered
    if (result != IGRAPH_SUCCESS && result != IGRAPH_STOP) {
        throw std::runtime_error("VF2 automorphism search failed with unexpected error");
    }
    
    // Return true if we found at least k automorphisms
    return ctx.current_count >= k || ctx.should_stop;
}


// Main interface: bypass + VF2 integration
bool AutomorphismCounter::hasAtLeastKAutomorphisms(const adjacency_matrix_t& matrix, int k) {
    // 1. Simple bypass (100% accurate, no false negatives)
    if (simpleBypass(matrix, k)) {
        return true;
    }
    
    // 2. VF2 with early termination (expensive fallback)
    return vf2EarlyTermination(matrix, k);
}

// Exact counting (for validation purposes)
long long AutomorphismCounter::countExactAutomorphisms(const adjacency_matrix_t& matrix) {
    // Check for unknown edges first
    bool hasUnknownEdges = false;
    for (vertex_t i = 0; i < vertices && !hasUnknownEdges; i++) {
        for (vertex_t j = 0; j < vertices && !hasUnknownEdges; j++) {
            if (matrix[i][j] == truth_value_unknown) {
                hasUnknownEdges = true;
            }
        }
    }
    
    if (hasUnknownEdges) {
        // For PDGs, use VF2 with exact counting (no early termination)
        igraph_t graph = convertToIgraph(matrix);
        
        // Set up VF2 context for exact counting
        VF2EarlyTerminationContext ctx;
        // Set target_count to a value that will never be reached,
        // ensuring all automorphisms are counted. Use a very large value instead of max to avoid conversion issues.
        ctx.target_count = 1000000000LL; // 1 billion should be more than enough for any practical case
        ctx.current_count = 0;
        ctx.matrix = &matrix; // Single matrix for automorphism
        ctx.should_stop = false; // Will remain false as target_count is effectively infinite
        
        // Use VF2 for automorphism counting
        igraph_error_t result = igraph_isomorphic_function_vf2(
            &graph, &graph,  // Same graph for automorphism
            nullptr, nullptr,  // No vertex colors
            nullptr, nullptr,  // No edge colors  
            nullptr, nullptr,  // No mapping output needed
            vf2_automorphism_found_callback,  // This callback will now count all
            nullptr,  // Default vertex compatibility (identity check for automorphism)
            vf2_edge_compatibility_callback,  // The fixed SMS edge compatibility callback
            &ctx  // Our context
        );
        
        destroyIgraph(graph);
        
        // IGRAPH_STOP is a valid return code when early termination is triggered (not expected here)
        if (result != IGRAPH_SUCCESS && result != IGRAPH_STOP) {
            throw std::runtime_error("VF2 exact automorphism search failed with unexpected error");
        }
        
        return ctx.current_count; // Returns the total count of automorphisms
    }
    
    // For fully defined graphs, use igraph's exact automorphism counting
    igraph_t graph = convertToIgraph(matrix);
    igraph_bliss_info_t info;
    memset(&info, 0, sizeof(info));
    
    igraph_error_t result = igraph_count_automorphisms(&graph, nullptr, IGRAPH_BLISS_FL, &info);
    destroyIgraph(graph);
    
    if (result != IGRAPH_SUCCESS) {
        if (info.group_size) {
            igraph_free(info.group_size);
        }
        throw std::runtime_error("igraph automorphism counting failed");
    }
    
    // Convert string result to long long
    long long count = 1;
    if (info.group_size) {
        count = strtoll(info.group_size, nullptr, 10);
        igraph_free(info.group_size);
    }
    
    return count;
}