/**
 * @file universal2.hpp
 * @author Markus Kirchweger
 *
 * This file contains the header for circuit based 2QBF solver based on a CEGAR approach.
 * It also implements gate hashing by ensuring that each gate is only created once.
 * The current implementation only allows And and Or gates.
 */

#include "../graphChecker.hpp"
#include "../cadical.hpp"
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
    // get a topological order of the gates considering polarity
    virtual void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited) = 0;

    // add encoding of the gate to cnf (not encoding subcircuits); dependent to which solver it is added
    virtual void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf) = 0;
    // add encoding of the gate to cnf (not encoding subcircuits) considering polarity
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
    void topologicalOrder(vector<GateOrLiteral *> &, std::unordered_set<GateOrLiteral *> &);
    void appendEncoding(solver_id, int &, cnf_t &);
    void topologicalOrderPolarity(polarity_t, vector<std::pair<polarity_t, GateOrLiteral *>> &, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &);
    void appendEncodingPolarity(polarity_t, solver_id, int &, cnf_t &);
    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &);
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

    bool isEqual(const GateOrLiteral *other) const;
    // still have to be consider in the topological order so that the variable is added to the solver, i.e., solver2gateVariable is set
    void topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited);
    void appendEncoding(solver_id solver, int &, cnf_t &);
    void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited);
    void appendEncodingPolarity(polarity_t p, solver_id solver, int &, cnf_t &);
    signedGate_t simplifyGate(const vector<truth_value_t> &assignment, std::unordered_map<GateOrLiteral *, signedGate_t> &);
    bool onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &);
};

class NAryGate : public GateOrLiteral
{
public:
    const vector<signedGate> &getInputs() { return inputs; }

protected:
    NAryGate(vector<signedGate> inputs);

    vector<signedGate> inputs;
    void topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited);
    void topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited);
    bool onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &visited);
};

class AndGate : public NAryGate
{
private:
    AndGate(vector<signedGate> inputs) : NAryGate(inputs) {}; // TODO compute hash (already done in Nary Gate but could adapt by gate type)
    static std::unordered_set<AndGate *, Hash, Equal> gates;  // contains all gates of this type which ever occured

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

    void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf);
    void appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf);
    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications);
    bool isEqual(const GateOrLiteral *other) const;
};

class OrGate : public NAryGate
{
private:
    OrGate(vector<signedGate> inputs) : NAryGate(inputs) {}; // TODO compute hash (already done in Nary Gate but could adapt by gate type)
    static std::unordered_set<OrGate *, Hash, Equal> gates;  // contains all gates of this type which ever occured

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

    void appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf);
    void appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf);
    signedGate_t simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications);
    bool isEqual(const GateOrLiteral *other) const;
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

// get the highest variable of the a qcir instance. (Not including the gates)
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
    /**
     * @brief Construct a new QCIRchecker object
     *
     * @param file File containing the qcir encoding
     * @param polarityHashing Indicates whether the polarity of the gates should be considered i.e., teitsin transformation only in a single direction
     */
    QCIRchecker(std::ifstream &file, bool polarityHashing = false);
    cnf_t getEncodingForGate(GateOrLiteral *g, solver_id s, int &nextFreeVariable, bool negated);
    void checkProperty(const adjacency_matrix_t &, const vector<int> &model, int &nextFreeVariable);
};