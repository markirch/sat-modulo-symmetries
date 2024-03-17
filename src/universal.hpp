#include "graphChecker.hpp"
#include "cadical.hpp"
#include <unordered_map>
#include <cassert>

class UniversalChecker : public ComplexFullyDefinedGraphChecker
{
    CaDiCaL::Solver *universalSolver;
    cnf_t forAllCNF;
    vector<int> forAllAsumptions; // variables which are assumptions

public:
    // file contains the encoding in dimacs format
    UniversalChecker(std::ifstream &file, vector<int> forAllAsumptions)
    {
        name = "UniversalChecker";
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

// #include "decisionTree.hpp"
#define functionLearning false // activate function learning similar to QFUN

class UniversalCheckerQCIR : public ComplexFullyDefinedGraphChecker
{
    CaDiCaL::Solver *universalSolver;
    cnf_t forAllCNF;
    vector<int> forAllAsumptions; // variables which are assumptions
    vector<int> gateVariables;
    const int indexOutputGate = 0; // TODO !!!!! ensure that this is the case

    bool hashGates = true;                                  // try to reuse some stuff over several calls
    vector<std::unordered_map<vector<bool>, int>> gate2int; // hash gates with specific inputs to avoid duplicates
    vector<vector<int>> gate2dependencies;                  // for each gate, get the universal variables it depends on.

    // vector<pair<vector<bool>, vector<bool>>> samples; // store the assignemnt of the existential variables and from the universal variables
    vector<int> universalVariables;
    // vector<int> existentialVariables;

    vector<bool> isGate;        // tells whether it is a gate or not
    vector<bool> isExistential; // tells whether it is a existential variable or not
    vector<int> gate2index;     // given a gate, get the corresponding index in the gateVariables vector

public:
    // file contains the encoding in dimacs format
    UniversalCheckerQCIR(std::ifstream &file, vector<int> forAllAsumptions)
    {
        name = "UniversalCheckerQCIR";
        int maxVar;
        file2cnf(file, forAllCNF, maxVar);
        this->forAllAsumptions = forAllAsumptions;

        for (auto clause : forAllCNF)
        {
            auto var = abs(clause[0]);
            if (find(gateVariables.begin(), gateVariables.end(), var) == gateVariables.end())
            {
                gateVariables.push_back(var);
            }
        }

        // printf("Gates:");
        // for (auto l : gateVariables)
        //     printf("%d ", l);
        // printf("\n");
        // fflush(stdout);

        universalSolver = new CaDiCaL::Solver();
        for (auto clause : forAllCNF)
        {
            for (auto lit : clause)
                universalSolver->add(lit);
            universalSolver->add(0);
        }
        // printf("Outputgate: %d\n", gateVariables[indexOutputGate]);

        if (hashGates)
        {
            int maxGate = *std::max_element(gateVariables.begin(), gateVariables.end());
            gate2int = vector<std::unordered_map<vector<bool>, int>>(maxGate + 1);
            gate2dependencies = vector<vector<int>>(maxGate + 1);
            gate2index = vector<int>(maxGate + 1, -1);

            assert(maxVar >= maxGate);
            isGate = vector<bool>(maxVar + 1, false);

            for (int i = 0; i < (int)gateVariables.size(); i++)
            {
                gate2index[gateVariables[i]] = i;
            }
            for (auto g : gateVariables)
                isGate[g] = true;
            isExistential = vector<bool>(maxVar + 1, false); // TODO tells whether it is an existential variable
            for (auto v : forAllAsumptions)
                isExistential[v] = true;

            int previousGate = -1;
            for (int i = 0; i < (int)forAllCNF.size(); i++)
            {
                auto clause = forAllCNF[i];
                int gate = abs(clause[0]);

                if (gate == gateVariables[indexOutputGate]) // TODO output gate is ignored for some special reasons. (comes first in encoding)
                    continue;

                if (gate == previousGate)
                    continue;
                previousGate = gate;

                for (int j = 1; j < (int)clause.size(); j++)
                {
                    int x = abs(clause[j]);
                    if (isGate[x])
                    {
                        auto d = gate2dependencies[x];
                        for (auto v : d)
                        {
                            if (std::find(gate2dependencies[gate].begin(), gate2dependencies[gate].end(), v) == gate2dependencies[gate].end())
                            {
                                gate2dependencies[gate].push_back(v); // add all the variables it is dependent from
                            }
                        }
                    }
                    else
                    {
                        if (!isExistential[x])
                        {
                            gate2dependencies[gate].push_back(x); // universal
                            if (std::find(universalVariables.begin(), universalVariables.end(), x) == universalVariables.end())
                            {
                                // If not present, add it
                                universalVariables.push_back(x);
                            }
                        }
                    }
                }

                // printf("Dependencies for gate %d: ", gate);
                // for (auto d : gate2dependencies[gate])
                //     printf("%d ", d);
                // printf("\n");
            }
        }
    }
    void checkProperty(const adjacency_matrix_t &matrix, const vector<int> &model, int &nextFreeVariable);
};
