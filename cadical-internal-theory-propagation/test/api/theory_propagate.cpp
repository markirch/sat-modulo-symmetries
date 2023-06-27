#include "../../src/cadical.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <iostream>
#include <cassert>
#include <cmath>
#include <set>
#include <deque>
#include <algorithm>

static int n = 13;

static CaDiCaL::Solver solver;

static int ph (int p, int h) {
    assert (0 <= p), assert (p < n + 1);
    assert (0 <= h), assert (h < n);
    return 1 + h * (n+1) + p;
}

static int pos (int row, int col) {
    assert (0 <= row), assert (row < n);
    assert (0 <= col), assert (col < n);
    return (n * row) + col + 1;
}

static void nQueensFormula(CaDiCaL::Solver & solver) {
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            solver.add (pos (i, j));
        solver.add(0);

        for (int j = 0; j < n; j++)
            solver.add (pos (j, i));
        solver.add(0);
    }

    for (int i = 0; i < n; i++)
        for (int j1 = 0; j1 < n-1; j1++)
            for (int j2 = j1 + 1; j2 < n; j2++) {
                // Exactly one queen per row
                solver.add (-pos (i, j1)),solver.add (-pos (i, j2)), solver.add (0);
                // Exactly one queen per column
                solver.add (-pos (j1, i)),solver.add (-pos (j2, i)), solver.add (0);
            }
  
    // at-most-one queen per Diagonal up and down
    for (int p = 0; p < n; p++) {
        int x, y, i, j;
        for (x = 0, y = p; x < n - 1 && y < n - 1; x++, y++)
            for (i = x + 1, j = y + 1; i < n && j < n; i++, j++) {
                solver.add (-pos (x, y)),solver.add (-pos (i, j)), solver.add (0);
        }
        for (x = p, y = 0; x < n - 1 && y < n - 1; x++, y++)
            for (i = x + 1, j = y + 1; i < n && j < n; i++, j++) {
                solver.add (-pos (x, y)),solver.add (-pos (i, j)), solver.add (0);
        }

        
        for (x = p, y = 0; x > 0 && y < n - 1; x--, y++)
            for (i = x - 1, j = y + 1; i >= 0 && j < n; i--, j++) {
                solver.add (-pos (x, y)),solver.add (-pos (i, j)), solver.add (0);
            }
        for (x = n - 1, y = p; x > 0 && y < n - 1; x--, y++)
            for (i = x - 1, j = y + 1; i >= 0 && j < n; i--, j++) {
                solver.add (-pos (x, y)),solver.add (-pos (i, j)), solver.add (0);
            }
    }
     
}

static void PHPformula (CaDiCaL::Solver & solver, bool satisfiable = false) {
    int p_limit = n + 1;
    if (satisfiable) p_limit = n;

    for (int h = 0; h < n; h++)
        for (int p1 = 0; p1 < p_limit; p1++)
            for (int p2 = p1 + 1; p2 < p_limit; p2++)
                solver.add (-ph (p1, h)), solver.add (-ph (p2, h)), solver.add (0);

    for (int p = 0; p < p_limit; p++) {
        for (int h = 0; h < n; h++)
            solver.add (ph (p, h));
        solver.add (0);
    }

    // Symmetry breaking clauses would be:
    // for (int p = 1; p < p_limit-1; p++) {
    //     for (int h = p + 1; h < n; h++) {
    //         solver.add (ph (p - 1, p)), solver.add (-ph (p - 1, h)), solver.add (0);
    //     }
    // }
}

static int pigeon_id (int lit) {
    int idx = abs (lit);
    return ((idx - 1) % (n + 1));
}

static int hole_id (int lit) {
    int idx = abs (lit);
    return (int)floor((idx - 1)/(n + 1));
}

static bool has_lit(const std::deque<std::vector<int>> & current_trail, int lit) {
    for (auto level_lits : current_trail) {
        if(std::find(level_lits.begin(), level_lits.end(), lit) != level_lits.end()) {
            return true;
        }
        
    }
    return false;
}

class SymmetryBreaker : CaDiCaL::ExternalPropagator {
    CaDiCaL::Solver * solver;
    std::vector<std::vector<int>> new_clauses;
    std::deque<std::vector<int>> current_trail;
    int count = 0;
public:
    SymmetryBreaker(CaDiCaL::Solver * s) : solver(s) { 
        solver->connect_external_propagator(this);
        for (int h = 0; h < n; h++)
            for (int p = 0; p < n + 1; p++)
            solver->add_observed_var(ph(p,h));

        // The root-level of the trail is always there
        current_trail.push_back(std::vector<int>());
    }

    ~SymmetryBreaker () {
        solver->disconnect_external_propagator (); 
    }

    void notify_assignment(int lit, bool is_fixed) {
        if (is_fixed) {
            current_trail.front().push_back(lit);
        } else {
            current_trail.back().push_back(lit);
        }

        int p = pigeon_id(lit); 
        int h = hole_id(lit); 
        int j = p + 1;
        assert(ph(p,h) == abs(lit));
        
        if (lit > 0) {
            if (p < 1 || p >= n - 1 || h < p + 1 || h >= n ) return;
            
            if (has_lit(current_trail,-ph(j-1,j))) {
                //std::cout << "c THEORY-LOG Conflict-1 is found over binary clause: " << ph(j-1,j) << " || " << -ph(j-1,h) << std::endl;
                std::vector<int> clause;
                clause.push_back(ph(j-1,j));
                clause.push_back(-ph(j-1,h));
                new_clauses.push_back(clause);
                return; 
            } else if (!has_lit(current_trail,ph(j-1,j))) {
                //std::cout << "c THEORY-LOG Propagation-1 is found over binary clause: " << ph(j-1,j) << " || " << -ph(j-1,h) << std::endl;
                std::vector<int> clause;               
                clause.push_back(ph(j-1,j));
                clause.push_back(-ph(j-1,h));        
                new_clauses.push_back(clause);
            }
        } else {
            if (p + 1 != h || p < 1 || p >= n - 1) return;
            assert (-lit == ph(p, p+1)); // ph(j-1, j) is assigned false --> -ph (j-1, k) must be non-falsified for all k > j
            
            for (int k = j + 1; k < n; k++) {
                if (has_lit(current_trail,ph(j-1,k))) {
                    //std::cout << "c THEORY-LOG Conflict-2 is found over binary clause: " << ph(j-1,j) << " || " << -ph(j-1,k) << std::endl;
                    std::vector<int> clause;
                    clause.push_back(ph(j-1,j));
                    clause.push_back(-ph(j-1,k));
                    new_clauses.push_back(clause);
                    return; 
                } else if (!has_lit(current_trail,-ph(j-1,k))) {
                    //std::cout << "c THEORY-LOG Propagation-2 is found over binary clause: " << ph(j-1,j) << " || " << -ph(j-1,k) << std::endl;
                    std::vector<int> clause;
                    clause.push_back(ph(j-1,j));
                    clause.push_back(-ph(j-1,k));
                    new_clauses.push_back(clause);
                }
            }
        }
    }

    
    void notify_new_decision_level () {
        current_trail.push_back(std::vector<int>()); 
    }
    
    void notify_backtrack (size_t new_level) {
        while (current_trail.size() > new_level + 1) {
            current_trail.pop_back();
        }
    }

    bool cb_check_found_model (const std::vector<int> & model) { 
        (void) model;    
        return true;
    }

    bool cb_has_external_clause () {
        return (!new_clauses.empty());
    }

    int cb_add_external_clause_lit () { 
        if (new_clauses.empty()) return 0;
        else {
            assert(!new_clauses.empty());
            size_t clause_idx = new_clauses.size() - 1;
            if (new_clauses[clause_idx].empty()) {
                new_clauses.pop_back();
                return 0;
            } 
            
            int lit = new_clauses[clause_idx].back();
            new_clauses[clause_idx].pop_back();
            return lit;
        }
    }

    int cb_decide () { return 0; }
    int cb_propagate () { return 0; }
    int cb_add_reason_clause_lit (int plit) {
        (void)plit;
        return 0;
    };
};


class SolutionEnumerator : CaDiCaL::ExternalPropagator {
    CaDiCaL::Solver * solver;
    std::vector<int> blocking_clause = {};
    size_t full_assignment_size = 0;
    
public:
    int sol_count = 0; //TODO check if the solutions are uniq

    SolutionEnumerator(CaDiCaL::Solver * s) : solver(s) { 
        solver->connect_external_propagator(this);
        for (int i = 0; i < n; i++) 
            for (int j = 0; j < n; j++)
                solver->add_observed_var (pos (i, j));
           
        is_lazy = true;
        full_assignment_size = solver->active ();
    }
    ~SolutionEnumerator () { solver->disconnect_external_propagator (); }

    void process_trail (const std::vector<int> & current_trail) {
        (void) current_trail;
        assert(!is_lazy || false); //if is_lazy is true, this function should never be called
        return;
    } 

    bool cb_check_found_model (const std::vector<int> & model) {
        assert(model.size() == full_assignment_size);
        sol_count += 1;

        blocking_clause.clear();

        // std::cout << "c New solution was found: ";
        for (const auto& lit: model) {
            if (lit > 0) {
                // std::cout << lit << " ";
                blocking_clause.push_back(-lit);
            }
        }
        // std::cout << std::endl;
        
        return false;
    }

    bool cb_has_external_clause () {
        return (!blocking_clause.empty());
    }

    int cb_add_external_clause_lit () { 
        if (blocking_clause.empty()) return 0;
        else {
            int lit = blocking_clause.back();
            blocking_clause.pop_back();
            return lit;
        }
    }



    void notify_assignment (int,bool) {};
    void notify_new_decision_level () {};
    void notify_backtrack (size_t) {};
    int cb_decide () { return 0; }
    int cb_propagate () { return 0; }
    int cb_add_reason_clause_lit (int) { return 0; };
};

void ph_test(bool satisfiable) {
    PHPformula(solver, satisfiable);
    int max_var = solver.active ();
    
    std::cout << "Solving PH-" << n << "." << std::endl;
    std::cout << "Nof vars: " << max_var << std::endl;

    SymmetryBreaker sb(&solver);

    int res = solver.solve ();
    std::cout << "Result: " << res << std::endl;
    assert (!satisfiable || res == 10);
    assert (satisfiable || res == 20);
    
    //solver.statistics ();
    solver.resources ();
}

void nQueens_test() {
    nQueensFormula(solver);
    int max_var = solver.active ();

    std::cout << "Enumerating solutions of the " << n << "-Queens problem." << std::endl;
    std::cout << "Nof vars: " << max_var << std::endl;

    SolutionEnumerator se(&solver);

    int res = solver.solve ();
    std::cout << "Result: " << res << std::endl;
    assert (res == 20);

    std::cout << "Number of enumerated solutions: " << se.sol_count << std::endl;

    solver.resources ();
}

int main () {
    solver.set("log",0);
    solver.set("chrono",0);
    solver.set("inprocessing",1);
    
    ph_test(false);
    //nQueens_test();
    //solver.statistics();
    return 0;
}