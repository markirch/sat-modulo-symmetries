#ifndef AUTOMORPHISM_COUNTER_HPP
#define AUTOMORPHISM_COUNTER_HPP

#include "../useful.h"
#include <vector>
#include <unordered_map>

class AutomorphismCounter {
private:
    int vertices;
    bool directed;
    
    // Nauty workspace members (allocated once for performance)
    void* g;          // graph* (avoiding nauty.h include in header)
    int* lab;
    int* ptn; 
    int* orbits;
    int m;            // SETWORDSNEEDED(vertices)
    size_t g_sz, lab_sz, ptn_sz, orbits_sz;  // For DYNFREE
    
    // Pure C++ helper methods
    int countTrueNeighbors(const adjacency_matrix_t& matrix, vertex_t v) const;
    int countZeroNeighborVertices(const adjacency_matrix_t& matrix) const;
    long long computeFactorial(int n) const;
    
    // Nauty FDG automorphism counting
    long long countFDGAutomorphismsNauty(const adjacency_matrix_t& matrix);
    
    // Nauty PDG threshold checking with enumeration
    bool nautyPDGThresholdCheck(const adjacency_matrix_t& matrix, int k);

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


#endif // AUTOMORPHISM_COUNTER_HPP