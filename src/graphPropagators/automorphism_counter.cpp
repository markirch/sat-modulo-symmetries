#include "automorphism_counter.hpp"
#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>
#include <vector>

extern "C" {
#include <nauty/nauty.h>
#include <nauty/naututil.h>
}

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


// Nauty PDG automorphism counting with enumeration (replacement for VF2)
bool AutomorphismCounter::nautyPDGThresholdCheck(const adjacency_matrix_t& matrix, int k) {
    // Find undefined edges
    vector<std::pair<int,int>> undefined_edges;
    if (directed) {
        // For directed graphs: check all (i,j) pairs where i != j
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < vertices; j++) {
                if (i != j && matrix[i][j] == truth_value_unknown) {
                    undefined_edges.push_back({i, j});
                }
            }
        }
    } else {
        // For undirected graphs: check upper triangle only
        for (int i = 0; i < vertices; i++) {
            for (int j = i+1; j < vertices; j++) {
                if (matrix[i][j] == truth_value_unknown) {
                    undefined_edges.push_back({i, j});
                }
            }
        }
    }
    
    int num_undefined = undefined_edges.size();
    
    // If no undefined edges, it's an FDG - use direct Nauty counting
    if (num_undefined == 0) {
        long long exact_count = countFDGAutomorphismsNauty(matrix);
        return exact_count >= k;
    }
    
    
    // For PDGs: enumerate all 2^k possible FDG assignments
    int total_assignments = 1 << num_undefined;
    
    for (int mask = 0; mask < total_assignments; mask++) {
        // Create FDG with this assignment
        auto fdg = matrix;
        for (int i = 0; i < num_undefined; i++) {
            int u = undefined_edges[i].first;
            int v = undefined_edges[i].second;
            truth_value_t value = (mask & (1 << i)) ? truth_value_true : truth_value_false;
            fdg[u][v] = value;
            if (!directed) {
                fdg[v][u] = value;  // Only set symmetric entry for undirected graphs
            }
        }
        
        // Count automorphisms for this FDG using Nauty
        long long auts = countFDGAutomorphismsNauty(fdg);
        
        // Early termination: if any FDG meets threshold, PDG can achieve it
        if (auts >= k) {
            return true;
        }
    }
    
    // No FDG assignment reached the threshold
    return false;
}

// Main interface: bypass + full Nauty integration  
bool AutomorphismCounter::hasAtLeastKAutomorphisms(const adjacency_matrix_t& matrix, int k) {
    // 1. Simple bypass (100% accurate, no false negatives)
    if (simpleBypass(matrix, k)) {
        return true;
    }
    
    // 2. Nauty for both FDGs and PDGs (complete VF2 replacement)
    return nautyPDGThresholdCheck(matrix, k);
}

// Nauty FDG automorphism counting (replacement for igraph/BLISS)
long long AutomorphismCounter::countFDGAutomorphismsNauty(const adjacency_matrix_t& matrix) {
    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    
    int m = SETWORDSNEEDED(vertices);
    
    DYNALLOC2(graph, g, g_sz, vertices, m, "malloc");
    DYNALLOC1(int, lab, lab_sz, vertices, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, vertices, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, vertices, "malloc");
    
    EMPTYGRAPH(g, m, vertices);
    
    // Configure options for directed vs undirected graphs
    options.getcanon = FALSE;
    options.defaultptn = TRUE;
    options.userautomproc = nullptr;
    options.userlevelproc = nullptr;
    options.digraph = directed ? TRUE : FALSE;
    
    // Add edges from FDG (only truth_value_true edges)
    if (directed) {
        // For directed graphs: use ADDONEARC and check all (i,j) pairs
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < vertices; j++) {
                if (i != j && matrix[i][j] == truth_value_true) {
                    ADDONEARC(g, i, j, m);
                }
            }
        }
    } else {
        // For undirected graphs: use ADDONEEDGE
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < vertices; j++) {
                if (matrix[i][j] == truth_value_true) {
                    ADDONEEDGE(g, i, j, m);
                }
            }
        }
    }
    
    densenauty(g, lab, ptn, orbits, &options, &stats, m, vertices, nullptr);
    
    long long total_automorphisms = stats.grpsize1;
    for (int i = 0; i < stats.grpsize2; i++) {
        total_automorphisms *= 10;
    }
    
    DYNFREE(g, g_sz);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    DYNFREE(orbits, orbits_sz);
    
    return total_automorphisms;
}

// Nauty exact counting for both FDGs and PDGs (complete VF2 replacement)
long long AutomorphismCounter::countExactAutomorphisms(const adjacency_matrix_t& matrix) {
    // Find undefined edges
    vector<std::pair<int,int>> undefined_edges;
    if (directed) {
        // For directed graphs: check all (i,j) pairs where i != j
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < vertices; j++) {
                if (i != j && matrix[i][j] == truth_value_unknown) {
                    undefined_edges.push_back({i, j});
                }
            }
        }
    } else {
        // For undirected graphs: check upper triangle only
        for (int i = 0; i < vertices; i++) {
            for (int j = i+1; j < vertices; j++) {
                if (matrix[i][j] == truth_value_unknown) {
                    undefined_edges.push_back({i, j});
                }
            }
        }
    }
    
    int num_undefined = undefined_edges.size();
    
    // If no undefined edges, it's an FDG - use direct Nauty counting
    if (num_undefined == 0) {
        return countFDGAutomorphismsNauty(matrix);
    }
    
    
    // For PDGs: find maximum automorphisms across all 2^k possible FDG assignments
    int total_assignments = 1 << num_undefined;
    long long max_automorphisms = 0;
    
    for (int mask = 0; mask < total_assignments; mask++) {
        // Create FDG with this assignment
        auto fdg = matrix;
        for (int i = 0; i < num_undefined; i++) {
            int u = undefined_edges[i].first;
            int v = undefined_edges[i].second;
            truth_value_t value = (mask & (1 << i)) ? truth_value_true : truth_value_false;
            fdg[u][v] = value;
            if (!directed) {
                fdg[v][u] = value;  // Only set symmetric entry for undirected graphs
            }
        }
        
        // Count automorphisms for this FDG using Nauty
        long long auts = countFDGAutomorphismsNauty(fdg);
        
        // Track maximum
        if (auts > max_automorphisms) {
            max_automorphisms = auts;
        }
    }
    
    return max_automorphisms;
}