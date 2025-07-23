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





// Nauty PDG automorphism counting with enumeration
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
    
    // 2. Nauty for both FDGs and PDGs
    return nautyPDGThresholdCheck(matrix, k);
}

// Nauty FDG automorphism counting
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

// Nauty exact counting for both FDGs and PDGs
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