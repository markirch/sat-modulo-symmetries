#include "graphChecker.hpp"
#include "cadical.hpp"

class UniversalChecker : public ComplexFullyDefinedGraphChecker
{
    CaDiCaL::Solver *universalSolver;
    cnf_t forAllCNF;
    vector<int> forAllAsumptions; // variables which are assumptions

public:
    // file contains the encoding in dimacs format
    UniversalChecker(std::ifstream &file, vector<int> forAllAsumptions)
    {
        int maxVar;
        file2cnf(file, forAllCNF, maxVar);
        this->forAllAsumptions = forAllAsumptions;


        universalSolver = new CaDiCaL::Solver();
        for (auto clause : forAllCNF)
        {
            for (auto lit : clause)
                universalSolver->add(lit);
            universalSolver->add(0);
        }
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable);
};