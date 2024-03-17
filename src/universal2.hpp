/**
 * @file universal2.hpp
 * @author Markus Kirchweger
 *
 * This file contains the header for circuit based 2QBF solver based on a CEGAR approach.
 * It also implements gate hashing by ensuring that each gate is only created once.
 * The current implementation only allows And and Or gates.
 */

#include "graphChecker.hpp"
#include "cadical.hpp"
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <boost/functional/hash.hpp>

typedef int solver_id; // used for identifying the solver. At the moment 1 is the existential and 2 the universal solver
typedef int polarity_t;

struct signedGate; // Forward declaration of signedGate_t

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};

class GateOrLiteral
{
protected:
    size_t hash; // must be computed in the constructor

public:
    virtual ~GateOrLiteral() {}

    std::unordered_map<solver_id, int> solver2gateVariable;                                              // map from solver to gate variable used for representing the gate
    std::unordered_map<std::pair<solver_id, polarity_t>, int, pair_hash> solverAndPolarity2gateVariable; // given solverid and polarity, return the gate variable if it exists. If the variable is present it also means that it is already encoded

    // get a topological order of the gates
    virtual void topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited) = 0;

    // add encoding to cnf; dependent to which solver it is added
    virtual void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf) = 0;

    // get a topological order of the gates considering polarity
    virtual void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited) = 0;
    virtual void appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf) = 0;

    // simplify gate based on given assignment
    virtual struct signedGate simplifyGate(const vector<truth_value_t> &assignment, std::unordered_map<GateOrLiteral *, struct signedGate> &simplifications) = 0;

    // check whether the gate is only dependent on the given variables otherwise return false
    virtual bool onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &visited) = 0;

    // check whether two gates are identical. Assumes that all previous generated gates can be identified by their pointers
    virtual bool isEqual(const GateOrLiteral *other) const = 0;

    struct Hash
    {
        size_t operator()(const GateOrLiteral *ptr) const
        {
            return ptr->hash;
        }
    };

    struct Equal
    {
        bool operator()(const GateOrLiteral *lhs, const GateOrLiteral *rhs) const
        {
            if (lhs->hash != rhs->hash)
                return false;
            return lhs->isEqual(rhs);
        }
    };

    size_t getHash() const { return hash; };
};

typedef struct signedGate
{
    int sign; // 1 if positive, -1 if negative
    GateOrLiteral *gate;
} signedGate_t;

// gate which is always true
class TrueGate : public GateOrLiteral
{
    void topologicalOrder(vector<GateOrLiteral *> &, std::unordered_set<GateOrLiteral *> &)
    {
        printf("Error: True gets should be removed before getting topological order\n");
        EXIT_UNWANTED_STATE
    }

    void appendEncoding(solver_id, int &, cnf_t &)
    {
        printf("Error: True gate should be removed before appending to encoding\n");
        EXIT_UNWANTED_STATE
    }

    void topologicalOrderPolarity(polarity_t, vector<std::pair<polarity_t, GateOrLiteral *>> &, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &)
    {
        printf("Error: True gets should be removed before getting topological order\n");
        EXIT_UNWANTED_STATE
    }

    void appendEncodingPolarity(polarity_t, solver_id, int &, cnf_t &)
    {
        printf("Error: True gate should be removed before appending to encoding\n");
        EXIT_UNWANTED_STATE
    }

    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &)
    {
        printf("Should never be called\n");
        EXIT_UNWANTED_STATE
    };

    bool onlyDependentOn(const std::unordered_set<int> &, std::unordered_set<GateOrLiteral *> &) { return true; }

public:
    static TrueGate *getInstance()
    {
        static TrueGate instance; // Static instance created once
        return &instance;
    }

    bool isEqual(const GateOrLiteral *other) const { return other == this; }; // only one object

private:
    // make it a singleton
    TrueGate() {}
};

class Variable : public GateOrLiteral
{
private:
    // constructor
    Variable(int variable) : variable(variable) { this->hash = std::hash<int>{}(variable); };
    static std::unordered_set<Variable *, Hash, Equal> gates; // contains all gates of this type which ever occured

public:
    const int variable;
    static Variable *getVariableGate(int variable)
    {
        Variable *v = new Variable(variable);
        auto it = gates.find(v);
        if (it != gates.end())
        {
            EXIT_UNWANTED_STATE // should not happen that a variable is created twice, must likely a bug in the code
                delete v;
            return *it;
        }
        Variable::gates.insert(v);
        return v;
    }

    bool isEqual(const GateOrLiteral *other) const
    {
        if (auto otherVariable = dynamic_cast<const Variable *>(other))
            return otherVariable->variable == variable;
        return false;
    }

    // still have to be consider in the topological order so that the variable is added to the solver, i.e., solver2gateVariable is set
    void topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited)
    {
        if (visited.find(this) != visited.end())
            return; // already visited
        visited.insert(this);
        order.push_back(this);
    }

    void appendEncoding(solver_id solver, int &, cnf_t &)
    {
        solver2gateVariable[solver] = variable; // TODO check that not negated
    }

    void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited)
    {
        if (visited.find({p, this}) != visited.end())
            return; // already visited
        visited.insert({p, this});
        order.push_back({p, this});
    }

    void appendEncodingPolarity(polarity_t p, solver_id solver, int &, cnf_t &)
    {
        solverAndPolarity2gateVariable[{solver, p}] = variable;
    }

    signedGate_t simplifyGate(const vector<truth_value_t> &assignment, std::unordered_map<GateOrLiteral *, signedGate_t> &)
    {
        if (assignment[variable] == truth_value_true)
            return {1, TrueGate::getInstance()};
        else if (assignment[variable] == truth_value_false)
            return {-1, TrueGate::getInstance()};
        else
            return {1, this};
    }

    bool onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &)
    {
        // no advantage to check whether it was already visited before
        return variables.find(variable) != variables.end();
    }
};

class NAryGate : public GateOrLiteral
{
public:
    const vector<signedGate> &getInputs() { return inputs; }

protected:
    vector<signedGate> inputs;
    void topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited)
    {
        if (visited.find(this) != visited.end())
            return; // already visited

        for (auto &input : inputs)
        {
            input.gate->topologicalOrder(order, visited);
        }
        visited.insert(this);
        order.push_back(this);
    }

    void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited)
    {
        if (visited.find({p, this}) != visited.end())
            return; // already visited

        for (auto &input : inputs)
        {
            input.gate->topologicalOrderPolarity(input.sign * p, order, visited);
        }
        visited.insert({p, this});
        order.push_back({p, this});
    }

    bool onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &visited)
    {
        // if already visited then skip to avoid blow up
        if (visited.find(this) != visited.end())
            return true;
        visited.insert(this);

        for (auto &input : inputs)
        {
            if (!input.gate->onlyDependentOn(variables, visited))
                return false;
        }
        return true;
    }

    NAryGate(vector<signedGate> inputs)
    {
        this->inputs = inputs;
        // compute hash
        hash = 0;
        for (auto &input : inputs)
        {
            size_t inputHash = input.gate->getHash();
            if (input.sign == -1)
                inputHash = ~inputHash;
            hash ^= inputHash + 0x9e3779b9 + (hash << 6) + (hash >> 2); // typical hash function; // TODO check clashes and improve if necessary
        }
    }
};

class AndGate : public NAryGate
{
private:
    AndGate(vector<signedGate> inputs) : NAryGate(inputs){}; // TODO compute hash (already done in Nary Gate but could adapt by gate type)
    static std::unordered_set<AndGate *, Hash, Equal> gates; // contains all gates of this type which ever occured

public:
    static AndGate *getAndGate(vector<signedGate> inputs)
    {
        AndGate tmpGate(inputs); // temporary to avoid touching and deleting heap to often
        auto it = gates.find(&tmpGate);
        if (it != gates.end())
        {
            return *it;
        }
        AndGate *g = new AndGate(inputs);
        AndGate::gates.insert(g);
        return g;
    }

    void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf)
    {
        if (solver2gateVariable.find(solver) != solver2gateVariable.end())
            return; // already present so also already encoded

        int variable = nextFreeVariable++;
        solver2gateVariable[solver] = variable;

        for (auto &input : inputs)
        {
            if (input.sign == 1)
                cnf.push_back({-variable, input.gate->solver2gateVariable[solver]});
            else
                cnf.push_back({-variable, -input.gate->solver2gateVariable[solver]});
        }
        clause_t clause = {variable};
        for (auto &input : inputs)
        {
            if (input.sign == 1)
                clause.push_back(-input.gate->solver2gateVariable[solver]);
            else
                clause.push_back(input.gate->solver2gateVariable[solver]);
        }
        cnf.push_back(clause);
    }

    void appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf)
    {
        if (solverAndPolarity2gateVariable.find({solver, p}) != solverAndPolarity2gateVariable.end())
            return; // already present so also already encoded

        if (p != 0 && solverAndPolarity2gateVariable.find({solver, 0}) != solverAndPolarity2gateVariable.end())
        {
            solverAndPolarity2gateVariable[{solver, p}] = solverAndPolarity2gateVariable[{solver, 0}];
            return;
        }

        bool pPosPresent = solverAndPolarity2gateVariable.find({solver, 1}) != solverAndPolarity2gateVariable.end();
        bool pNegPresent = solverAndPolarity2gateVariable.find({solver, -1}) != solverAndPolarity2gateVariable.end();

        int variable; // reuse variable if already one is present, no matter the polarity
        if (pPosPresent)
            variable = solverAndPolarity2gateVariable[{solver, 1}];
        else if (pNegPresent)
            variable = solverAndPolarity2gateVariable[{solver, -1}];
        else
            variable = nextFreeVariable++;

        solverAndPolarity2gateVariable[{solver, p}] = variable;

        if (p == 1 || (p == 0 && !pPosPresent))
        {
            for (auto &input : inputs)
            {
                if (input.sign == 1)
                    cnf.push_back({-variable, input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]});
                else
                    cnf.push_back({-variable, -input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]});
            }
        }

        if (p == -1 || (p == 0 && !pNegPresent))
        {
            clause_t clause = {variable};
            for (auto &input : inputs)
            {
                if (input.sign == 1)
                    clause.push_back(-input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]);
                else
                    clause.push_back(input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]);
            }
            cnf.push_back(clause);
        }
    }

    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications)
    {
        vector<signedGate_t> simplifiedInputs;

        bool changeInInputs = false;
        for (auto &input : inputs)
        {
            signedGate_t simplifiedInput = simplifications[input.gate]; // must already been simplified
            if (simplifiedInput.sign != 1 || simplifiedInput.gate != input.gate)
                changeInInputs = true;

            // change sign depending on input
            simplifiedInput = {simplifiedInput.sign * input.sign, simplifiedInput.gate};

            if (simplifiedInput.gate == TrueGate::getInstance())
            {
                if (simplifiedInput.sign == 1)
                    continue; // input gate is always true
                else
                    return {-1, TrueGate::getInstance()}; // gate is always false, i.e., overall gate is always false
            }

            // base case
            simplifiedInputs.push_back(simplifiedInput);
        }

        if (simplifiedInputs.size() == 0)
            return {1, TrueGate::getInstance()};

        if (simplifiedInputs.size() == 1)
            return simplifiedInputs[0];

        if (changeInInputs)
        {
            return {1, getAndGate(simplifiedInputs)}; // TODO check whether already hashed before creating a new one.
        }
        return {1, this};
    };

    bool isEqual(const GateOrLiteral *other) const // TODO eventually replace in Nary gate and use type id
    {
        if (auto otherAndGate = dynamic_cast<const AndGate *>(other))
        {

            for (auto i = 0; i < (int)inputs.size(); i++) // all inputs must be equal and have the same sign
            {

                if (inputs[i].sign != otherAndGate->inputs[i].sign ||
                    inputs[i].gate != otherAndGate->inputs[i].gate) // TODO check again but can use pointer because of using singeltons
                {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
};

/**
 * @brief
 *
 */
class OrGate : public NAryGate
{
private:
    OrGate(vector<signedGate> inputs) : NAryGate(inputs){}; // TODO compute hash (already done in Nary Gate but could adapt by gate type)
    static std::unordered_set<OrGate *, Hash, Equal> gates; // contains all gates of this type which ever occured

public:
    static OrGate *getOrGate(vector<signedGate> inputs)
    {
        OrGate tmpGate(inputs); // temporary to avoid touching and deleting heap to often
        auto it = gates.find(&tmpGate);
        if (it != gates.end())
        {
            return *it;
        }
        OrGate *g = new OrGate(inputs);
        OrGate::gates.insert(g);
        return g;
    }

    void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf)
    {
        if (solver2gateVariable.find(solver) != solver2gateVariable.end())
            return; // already present so also already encoded

        int variable = nextFreeVariable++;
        solver2gateVariable[solver] = variable;

        for (auto &input : inputs)
        {
            if (input.sign == 1)
                cnf.push_back({variable, -input.gate->solver2gateVariable[solver]});
            else
                cnf.push_back({variable, input.gate->solver2gateVariable[solver]});
        }
        clause_t clause = {-variable};
        for (auto &input : inputs)
        {
            if (input.sign == 1)
                clause.push_back(input.gate->solver2gateVariable[solver]);
            else
                clause.push_back(-input.gate->solver2gateVariable[solver]);
        }
        cnf.push_back(clause);
    }

    void appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf)
    {
        if (solverAndPolarity2gateVariable.find({solver, p}) != solverAndPolarity2gateVariable.end())
            return; // already present so also already encoded

        if (p != 0 && solverAndPolarity2gateVariable.find({solver, 0}) != solverAndPolarity2gateVariable.end())
        {
            solverAndPolarity2gateVariable[{solver, p}] = solverAndPolarity2gateVariable[{solver, 0}];
            return;
        }

        bool pPosPresent = solverAndPolarity2gateVariable.find({solver, 1}) != solverAndPolarity2gateVariable.end();
        bool pNegPresent = solverAndPolarity2gateVariable.find({solver, -1}) != solverAndPolarity2gateVariable.end();

        int variable; // reuse variable if already one is present, no matter the polarity
        if (pPosPresent)
            variable = solverAndPolarity2gateVariable[{solver, 1}];
        else if (pNegPresent)
            variable = solverAndPolarity2gateVariable[{solver, -1}];
        else
            variable = nextFreeVariable++;

        solverAndPolarity2gateVariable[{solver, p}] = variable;

        if (p == 1 || (p == 0 && !pPosPresent))
        {
            clause_t clause = {-variable};
            for (auto &input : inputs)
            {
                if (input.sign == 1)
                    clause.push_back(input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]);
                else
                    clause.push_back(-input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]);
            }
            cnf.push_back(clause);
        }

        if (p == -1 || (p == 0 && !pNegPresent))
        {
            for (auto &input : inputs)
            {
                if (input.sign == 1)
                    cnf.push_back({variable, -input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]});
                else
                    cnf.push_back({variable, input.gate->solverAndPolarity2gateVariable[{solver, input.sign * p}]});
            }
        }
    }

    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications)
    {
        vector<signedGate_t> simplifiedInputs;

        bool changeInInputs = false;
        for (auto &input : inputs)
        {
            signedGate_t simplifiedInput = simplifications[input.gate]; // must already been simplified
            if (simplifiedInput.sign != 1 || simplifiedInput.gate != input.gate)
                changeInInputs = true;

            // change sign depending on input
            simplifiedInput = {simplifiedInput.sign * input.sign, simplifiedInput.gate};

            if (simplifiedInput.gate == TrueGate::getInstance())
            {
                if (simplifiedInput.sign == 1)
                    return {1, TrueGate::getInstance()}; // or-gate is always true
                else
                    continue; // input always false
            }

            // base case
            simplifiedInputs.push_back(simplifiedInput);
        }

        if (simplifiedInputs.size() == 0)
            return {-1, TrueGate::getInstance()}; // Empty or-gate is always false by definition

        if (simplifiedInputs.size() == 1)
            return simplifiedInputs[0];

        if (changeInInputs)
        {
            return {1, getOrGate(simplifiedInputs)}; // TODO check whether already hashed before creating a new one.
        }

        return {1, this};
    };

    bool isEqual(const GateOrLiteral *other) const
    {
        if (auto otherOrGate = dynamic_cast<const OrGate *>(other))
        {
            for (auto i = 0; i < (int)inputs.size(); i++) // all inputs must be equal and have the same sign
            {

                if (inputs[i].sign != otherOrGate->inputs[i].sign ||
                    inputs[i].gate != otherOrGate->inputs[i].gate) // TODO check again but can use pointer because of using singeltons
                {
                    return false;
                }
            }
            return true;
        }
        return false;
    }
};

typedef struct
{
    GateOrLiteral *outputGate;
    vector<int> freeVariables;
    vector<std::pair<int, vector<int>>> existentialVariables; // first element gives the quantifier depth, second the variables
    vector<std::pair<int, vector<int>>> universalVariables;   // first element gives the quantifier depth, second the variables
} qcir_t;

qcir_t parseQCIRFormula(std::ifstream &input);
void checkInput(qcir_t &instance);

// assuming that the first gate is an AndGate, we remove all input gates which only depend on existential variables.
GateOrLiteral *removeExistentialPart(GateOrLiteral *outputGate, vector<int> &existentialAndFreeVariables);

inline int getHighestVariable(qcir_t &instance)
{
    int highestVariable = 0;
    for (auto &universalPart : instance.universalVariables)
        for (auto &variable : universalPart.second)
            highestVariable = std::max(highestVariable, variable);
    for (auto &existentialPart : instance.existentialVariables)
        for (auto &variable : existentialPart.second)
            highestVariable = std::max(highestVariable, variable);
    for (auto &variable : instance.freeVariables)
        highestVariable = std::max(highestVariable, variable);
    return highestVariable;
}

// simplify gate based on given assignment
signedGate_t simplifyGate(GateOrLiteral *g, const vector<truth_value_t> &assignment);
// create cnf encoding of gate
cnf_t getEncodingForGate(GateOrLiteral *g, solver_id s, int &nextFreeVariable, bool negated);

inline void printCNF(cnf_t &cnf)
{
    for (auto clause : cnf)
    {
        for (auto literal : clause)
        {
            printf("%d ", literal);
        }
        printf("\n");
    }
}

class QCIRchecker : public ComplexFullyDefinedGraphChecker
{
    CaDiCaL::Solver *universalSolver;
    qcir_t instance;

    vector<int> existentialVariables; // for the sake of simplicity merge existential and free variables
    vector<int> universalVariables;

    bool existentialPartAdded = false;
    bool polarityHashing = false;

public:
    int highestVariableInInstance;

    // file contains the encoding in dimacs format
    QCIRchecker(std::ifstream &file, bool polarityHashing = false) : polarityHashing(polarityHashing)
    {
        name = "UniversalChecker";
        universalSolver = new CaDiCaL::Solver(); // id 2
        instance = parseQCIRFormula(file);

        checkInput(instance);

        // set existential variables based on instance
        if (instance.existentialVariables.size() > 0)
            existentialVariables = instance.existentialVariables[0].second;
        existentialVariables.insert(existentialVariables.end(), instance.freeVariables.begin(), instance.freeVariables.end());
        // set universal variables based on instance
        if (instance.universalVariables.size() > 0)
            universalVariables = instance.universalVariables[0].second;

        highestVariableInInstance = getHighestVariable(instance);
    }

    cnf_t getEncodingForGate(GateOrLiteral *g, solver_id s, int &nextFreeVariable, bool negated);

    void checkProperty(const adjacency_matrix_t &, const vector<int> &model, int &nextFreeVariable)
    {
        if (!existentialPartAdded)
        {
            existentialPartAdded = true;

            solver_id s = 1;

            // simplify with arbitrary universal assignment, i.e., setting all universal variables to false
            vector<truth_value_t> assignment(highestVariableInInstance + 1, truth_value_unknown);
            for (auto v : universalVariables)
                assignment[v] = truth_value_false;
            signedGate_t g = simplifyGate(instance.outputGate, assignment);

            cnf_t cnf = getEncodingForGate(g.gate, s, nextFreeVariable, g.sign == -1); // add gate; negate if gate has negative sign
            cnf.push_back({highestVariableInInstance, -highestVariableInInstance}); // just to ensure that all variables are set in the solver

            printf("Number of clauses added after first iteration of 2QBF solver: %ld\n", cnf.size());

            // remove existential part
            instance.outputGate = removeExistentialPart(instance.outputGate, existentialVariables);

            // load formula into universal solver
            int nextFreeVariableUniversalSolver = highestVariableInInstance + 1;
            s = 2; // so it is added to the universal solver
            // TODO eventually simplify before
            cnf_t cnfUniversal = getEncodingForGate(instance.outputGate, s, nextFreeVariableUniversalSolver, true); // negate output
            cnfUniversal.push_back({highestVariableInInstance, -highestVariableInInstance});                        // just to ensure that all variables are set in the solver

            // printf("Universel formula:\n");
            // printCNF(cnfUniversal);

            // load negated formula into universal solver
            for (auto clause : cnfUniversal)
            {
                for (auto literal : clause)
                {
                    universalSolver->add(literal);
                }
                universalSolver->add(0);
            }

            // printf("CNF for existential part:\n");
            // printCNF(cnf);

            throw cnf; // throw formula adding existental part
        }

        // printf("Highest variable in instance: %d, %ld\n", highestVariableInInstance, model.size());

        // check based on universal assignment given existential and free as assumptions
        for (auto v : existentialVariables)
        {
            universalSolver->assume(model[v - 1]);
        }
        // printf("Existential model: ");
        // for (auto v : existentialVariables)
        // {
        //     printf("(%d,%d) ", v, model[v - 1]);
        // }
        // printf("\n");

        auto res = universalSolver->solve();
        if (res != 10)
        {
            // printf("no counter model found, i.e., solution\n");
            return;
        }
        vector<truth_value_t> assignment(highestVariableInInstance + 1, truth_value_unknown);
        for (auto v : universalVariables)
            assignment[v] = universalSolver->val(v) > 0 ? truth_value_true : truth_value_false;

        // checkWhetherValidSpanningTree(matrix, universalVariables, universalSolver);
        // printf("Counter model: ");
        // for (auto v : universalVariables)
        // {
        //     printf("(%d,%d) ", v, universalSolver->val(v));
        // }
        // printf("\n");

        // TODO might also add root level fixed variables, but better to simplify the formula before

        signedGate_t g = simplifyGate(instance.outputGate, assignment);

        if (g.gate == TrueGate::getInstance() && g.sign == 1)
        {
            printf("Counter makes formula true\n");
            EXIT_UNWANTED_STATE
        }

        // add to encoding
        solver_id s = 1;
        cnf_t cnf = getEncodingForGate(g.gate, s, nextFreeVariable, g.sign == -1); // add gate; negate if gate has negative sign

        // printf("Clause from countermodel:");
        // for (clause_t clause : cnf)
        // {
        //     for (auto literal : clause)
        //     {
        //         printf("%d ", literal);
        //     }
        //     printf("\n");
        // }

        bool learnPertubations = false;
        if (learnPertubations)
        {
            printf("Number of universal variables: %ld\n", universalVariables.size());
            for (auto v : universalVariables)
            {
                // if (c++ > 10)
                //     break;
                assignment[v] = assignment[v] == truth_value_true ? truth_value_false : truth_value_true;

                signedGate_t g = simplifyGate(instance.outputGate, assignment);
                cnf_t pertubation = getEncodingForGate(g.gate, s, nextFreeVariable, g.sign == -1);
                cnf.insert(cnf.end(), pertubation.begin(), pertubation.end());

                assignment[v] = assignment[v] == truth_value_true ? truth_value_false : truth_value_true;
            }
        }

        throw cnf;
    }
};

// inline void checkWhetherValidSpanningTree(const adjacency_matrix_t &matrix, const vector<int> universalVariables,  CaDiCaL::Solver *universalSolver)
// {
//     // extract spanning tree from universal variables
//         int n = matrix.size();
//         adjacency_matrix_t spanningTree(n, vector<truth_value_t>(n, truth_value_unknown));
//         int c = 0;
//         for (int i = 0; i < n; i++)
//         {
//             for (int j = i + 1; j < n; j++)
//             {
//                 if (matrix[i][j] == 1)
//                 {
//                     spanningTree[i][j] = spanningTree[j][i] = universalSolver->val(universalVariables[c]) > 0 ? truth_value_true : truth_value_false;
//                 }
//                 c++;
//             }
//         }
//         int edgeCount = 0;
//         for (int i = 0; i < n; i++)
//             for (int j = i + 1; j < n; j++)
//             {
//                 if (spanningTree[i][j] == truth_value_true && matrix[i][j] == truth_value_false)
//                 {
//                     printf("Error: Spanning tree not a subgraph\n");
//                     EXIT_UNWANTED_STATE
//                 }
//                 if (spanningTree[i][j] == truth_value_true)
//                     edgeCount++;
//             }

//         if (edgeCount != n - 1)
//         {
//             printf("Error: Spanning tree has not enough edges\n");
//             EXIT_UNWANTED_STATE
//         }

//         // check that spanning tree is connected
//         vector<bool> visited(n, false);
//         visited[0] = true;
//         for (int round = 0; round < n * n; round++)
//             for (int i = 0; i < n; i++)
//             {
//                 if (visited[i])
//                 {
//                     for (int j = 0; j < n; j++)
//                     {
//                         if (spanningTree[i][j] == truth_value_true)
//                         {
//                             visited[j] = true;
//                         }
//                     }
//                 }
//             }
//         for (int i = 0; i < n; i++)
//         {
//             if (!visited[i])
//             {
//                 printf("Error: Spanning tree not connected\n");
//                 // print edges
//                 for (int i = 0; i < n; i++)
//                     for (int j = i + 1; j < n; j++)
//                     {
//                         if (spanningTree[i][j] == truth_value_true)
//                             printf("(%d,%d) ", i, j);
//                     }

//                 EXIT_UNWANTED_STATE
//             }
//         }

//         // check graph without spanning tree
//         vector<int> remainingDegrees(n, 0);
//         for (int i = 0; i < n; i++)
//             for (int j = i + 1; j < n; j++)
//             {
//                 if (matrix[i][j] == truth_value_true && spanningTree[i][j] == truth_value_false)
//                 {
//                     remainingDegrees[i]++;
//                     remainingDegrees[j]++;
//                 }
//             }
//         for (int i = 0; i < n; i++)
//             for (int j = 0; j < n; j++)
//             {
//                 if (matrix[i][j] == truth_value_true && spanningTree[i][j] == truth_value_false)
//                 {
//                     if (remainingDegrees[i] == 1 && remainingDegrees[j] == 2)
//                     {
//                         printf("Error: no valid spanning tree\n");
//                         // print spanning tree and matrix
//                         printf("Spanning tree:\n");
//                         for (int i = 0; i < n; i++)
//                         {
//                             for (int j = i + 1; j < n; j++)
//                             {
//                                 if (spanningTree[i][j] == truth_value_true)
//                                     printf("(%d,%d) ", i, j);
//                             }
//                         }
//                         printf("\n");

//                         printf("Matrix:\n");
//                         for (int i = 0; i < n; i++)
//                         {
//                             for (int j = i + 1; j < n; j++)
//                             {
//                                 if (matrix[i][j] == truth_value_true)
//                                     printf("(%d,%d) ", i, j);
//                             }
//                         }
//                         printf("\n");
//                         EXIT_UNWANTED_STATE
//                     }
//                 }
//             }

// }
