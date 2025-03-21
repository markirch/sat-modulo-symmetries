#include "universal2.hpp"

// Declare all static variables used in classes
std::unordered_set<Variable *, GateOrLiteral::Hash, GateOrLiteral::Equal> Variable::gates;
std::unordered_set<AndGate *, GateOrLiteral::Hash, GateOrLiteral::Equal> AndGate::gates;
std::unordered_set<OrGate *, GateOrLiteral::Hash, GateOrLiteral::Equal> OrGate::gates;

// simplify gate based on given assignment
signedGate_t simplifyGate(GateOrLiteral *g, const vector<truth_value_t> &assignment)
{
    vector<GateOrLiteral *> topologicalOrder;
    std::unordered_set<GateOrLiteral *> visited;
    g->topologicalOrder(topologicalOrder, visited);

    std::unordered_map<GateOrLiteral *, signedGate_t> simplifiedGates;
    // by using the topological order, we can ensure that all gates on which a gate is depending are already simplified.

    // for (auto gate : topologicalOrder)
    // {
    //     PRINT_CURRENT_LINE
    //     gate->isEqual(gate);
    //     PRINT_CURRENT_LINE
    // }

    for (auto gate : topologicalOrder)
    {
        // print type of gate
        // printf("Gate type: %s\n", typeid(*gate).name());
        // fflush(stdout);
        signedGate_t s = gate->simplifyGate(assignment, simplifiedGates);
        simplifiedGates[gate] = s;
    }

    return simplifiedGates[g];
}

// assuming that the first gate is an AndGate, we remove all input gates which only depend on existential variables.
GateOrLiteral *removeExistentialPart(GateOrLiteral *outputGate, vector<int> &existentialAndFreeVariables)
{
    if (auto andGate = dynamic_cast<AndGate *>(outputGate))
    {
        // transform vector into set for faster lookup
        std::unordered_set<int> existentialAndFreeVariablesSet(existentialAndFreeVariables.begin(), existentialAndFreeVariables.end());
        int numRemoved = 0;
        vector<signedGate_t> newInputGates;
        for (auto input : andGate->getInputs())
        {
            std::unordered_set<GateOrLiteral *> alreadyVisited;
            if (!input.gate->onlyDependentOn(existentialAndFreeVariablesSet, alreadyVisited))
            {
                newInputGates.push_back(input);
            }
            else
            {
                numRemoved++;
            }
        }

        if (numRemoved == 0)
            return outputGate;

        printf("Removed %d gates from AndGate\n", numRemoved);
        return AndGate::getAndGate(newInputGates);
    }
    else
    {
        return outputGate;
    }
}

// ------------------TrueGate------------------

void TrueGate::topologicalOrder(vector<GateOrLiteral *> &, std::unordered_set<GateOrLiteral *> &)
{
    printf("Error: True gets should be removed before getting topological order\n");
    EXIT_UNWANTED_STATE
}

void TrueGate::appendEncoding(solver_id, int &, cnf_t &)
{
    printf("Error: True gate should be removed before appending to encoding\n");
    EXIT_UNWANTED_STATE
}

void TrueGate::topologicalOrderPolarity(polarity_t, vector<std::pair<polarity_t, GateOrLiteral *>> &, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &)
{
    printf("Error: True gets should be removed before getting topological order\n");
    EXIT_UNWANTED_STATE
}

void TrueGate::appendEncodingPolarity(polarity_t, solver_id, int &, cnf_t &)
{
    printf("Error: True gate should be removed before appending to encoding\n");
    EXIT_UNWANTED_STATE
}

signedGate_t TrueGate::simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &)
{
    printf("Should never be called\n");
    EXIT_UNWANTED_STATE
}

// ------------------Variable------------------

bool Variable::isEqual(const GateOrLiteral *other) const
{
    if (auto otherVariable = dynamic_cast<const Variable *>(other))
        return otherVariable->variable == variable;
    return false;
}

// still have to be consider in the topological order so that the variable is added to the solver, i.e., solver2gateVariable is set
void Variable::topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited)
{
    if (visited.find(this) != visited.end())
        return; // already visited
    visited.insert(this);
    order.push_back(this);
}

void Variable::appendEncoding(solver_id solver, int &, cnf_t &)
{
    solver2gateVariable[solver] = variable; // TODO check that not negated
}

void Variable::topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited)
{
    if (visited.find({p, this}) != visited.end())
        return; // already visited
    visited.insert({p, this});
    order.push_back({p, this});
}

void Variable::appendEncodingPolarity(polarity_t p, solver_id solver, int &, cnf_t &)
{
    solverAndPolarity2gateVariable[{solver, p}] = variable;
}

signedGate_t Variable::simplifyGate(const vector<truth_value_t> &assignment, std::unordered_map<GateOrLiteral *, signedGate_t> &)
{
    if (assignment[variable] == truth_value_true)
        return {1, TrueGate::getInstance()};
    else if (assignment[variable] == truth_value_false)
        return {-1, TrueGate::getInstance()};
    else
        return {1, this};
}

bool Variable::onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &)
{
    // no advantage to check whether it was already visited before
    return variables.find(variable) != variables.end();
}

// ------------------NAryGate------------------

NAryGate::NAryGate(vector<signedGate> inputs)
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

void NAryGate::topologicalOrder(vector<GateOrLiteral *> &order, std::unordered_set<GateOrLiteral *> &visited)
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

void NAryGate::topologicalOrderPolarity(polarity_t p, vector<std::pair<polarity_t, GateOrLiteral *>> &order, std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> &visited)
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

bool NAryGate::onlyDependentOn(const std::unordered_set<int> &variables, std::unordered_set<GateOrLiteral *> &visited)
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

// ------------------AndGate------------------

void AndGate::appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf)
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

void AndGate::appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf)
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

signedGate_t AndGate::simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications)
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

bool AndGate::isEqual(const GateOrLiteral *other) const // TODO eventually replace in Nary gate and use type id
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

// ------------------OrGate------------------

void OrGate::appendEncoding(solver_id solver, int &nextFreeVariable, cnf_t &cnf)
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

void OrGate::appendEncodingPolarity(polarity_t p, solver_id solver, int &nextFreeVariable, cnf_t &cnf)
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

signedGate_t OrGate::simplifyGate(const vector<truth_value_t> &, std::unordered_map<GateOrLiteral *, signedGate_t> &simplifications)
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

bool OrGate::isEqual(const GateOrLiteral *other) const
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

// ----------- QCIRchecker ------------

QCIRchecker::QCIRchecker(std::ifstream &file, bool polarityHashing) : polarityHashing(polarityHashing)
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

cnf_t QCIRchecker::getEncodingForGate(GateOrLiteral *g, solver_id s, int &nextFreeVariable, bool negated)
{
    if (g == TrueGate::getInstance())
    {
        if (negated)
            return cnf_t({{-1}, {1}});
        else
            return cnf_t({{1, -1}});
    }

    if (polarityHashing)
    {
        vector<std::pair<polarity_t, GateOrLiteral *>> topologicalOrder;
        std::unordered_set<std::pair<polarity_t, GateOrLiteral *>, pair_hash> visited;
        g->topologicalOrderPolarity(negated ? -1 : 1, topologicalOrder, visited);

        cnf_t cnf;
        for (auto polarityGate : topologicalOrder)
        {
            auto gate = polarityGate.second;
            gate->appendEncodingPolarity(polarityGate.first, s, nextFreeVariable, cnf);
        }
        //  ensure gate itself is added
        int gateVariable = negated ? g->solverAndPolarity2gateVariable[{s, -1}] : g->solverAndPolarity2gateVariable[{s, 1}];
        cnf.push_back({negated ? -gateVariable : gateVariable});
        return cnf;
    }
    else
    {

        vector<GateOrLiteral *> topologicalOrder;
        std::unordered_set<GateOrLiteral *> visited;
        g->topologicalOrder(topologicalOrder, visited);

        cnf_t cnf;
        for (auto gate : topologicalOrder)
        {
            gate->appendEncoding(s, nextFreeVariable, cnf);
        }
        //  ensure gate itself is added
        int gateVariable = g->solver2gateVariable[s];
        cnf.push_back({negated ? -gateVariable : gateVariable});
        return cnf;
    }
}

void QCIRchecker::checkProperty(const adjacency_matrix_t &, const vector<int> &model, int &nextFreeVariable)
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
        cnf.push_back({highestVariableInInstance, -highestVariableInInstance});    // just to ensure that all variables are set in the solver

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
        assert(abs(model[v]) == v);
        universalSolver->assume(model[v]);
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
};
