#include "universal2.hpp"

// Declare all static variables used in classes
std::unordered_set<Variable *, GateOrLiteral::Hash, GateOrLiteral::Equal> Variable::gates;
std::unordered_set<AndGate *, GateOrLiteral::Hash, GateOrLiteral::Equal> AndGate::gates;
std::unordered_set<OrGate *, GateOrLiteral::Hash, GateOrLiteral::Equal> OrGate::gates;

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

/**

// get purely existential part from gate given the existential (and free) variables  from the first level.
GateOrLiteral *getExistentialPart(GateOrLiteral *outputGate, vector<int> &existentialVariables)
{
}

*/

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
