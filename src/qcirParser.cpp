#include "universal2.hpp"

// analog to python's strip
std::string trim(const std::string &str)
{
    std::string s = str;
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
                                    { return !std::isspace(ch); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
                         { return !std::isspace(ch); })
                .base(),
            s.end());
    return s;
}

// comma seperated integers to vector of integers
inline vector<int> line2ints(string line)
{
    vector<int> res;
    std::istringstream iss(line);
    clause_t clause;

    string lit;
    while (std::getline(iss, lit, ','))
    {
        int l = stoi(trim(lit));
        res.push_back(l);
    }
    return res;
}

// return true if line starts with prefix and parse the rest of the line
inline bool parseLinePrefix(const std::string &line, const std::string &prefix, std::vector<int> &res)
{
    if (line.substr(0, prefix.size()) == prefix)
    {
        res.clear(); // just to be sure
        if (line[prefix.size()] != '(' || line[line.size() - 1] != ')')
        {
            printf("Error: Missing opening or closing bracket for line %s\n", line.c_str());
            EXIT_UNWANTED_STATE
        }
        res = line2ints(line.substr(prefix.size() + 1, line.size() - (prefix.size() + 1 + 1))); // ignore brackets and quantifier
        return true;
    }
    return false;
}

// parse formula in QCIR format and extract circuit
qcir_t parseQCIRFormula(std::ifstream &input)
{
    qcir_t res;

    // get first line
    std::string line;
    std::getline(input, line);
    if (line != "#QCIR-G14")
    {
        printf("Error: File is not in QCIR-G14 format\n");
        EXIT_UNWANTED_STATE
    }

    // get quantification of the variables
    vector<string> quantifierPrefixes = {"free", "exists", "forall"};
    int quantifierDepth = 0;
    while (std::getline(input, line))
    {
        std::vector<int> variables;
        if (parseLinePrefix(line, "free", variables))
        {
            res.freeVariables = variables;
        }
        else if (parseLinePrefix(line, "exists", variables))
        {
            res.existentialVariables.push_back({quantifierDepth, variables});
            quantifierDepth++;
        }
        else if (parseLinePrefix(line, "forall", variables))
        {
            res.universalVariables.push_back({quantifierDepth, variables});
            quantifierDepth++;
        }
        else
        {
            break;
        }
    }

    int outputGateName;
    vector<int> inputs;
    parseLinePrefix(line, "output", inputs);
    if (inputs.size() != 1)
    {
        printf("Error: Only one output gate can be specified\n");
        EXIT_UNWANTED_STATE
    }
    outputGateName = inputs[0];

    // map from gateName to gate
    std::unordered_map<int, GateOrLiteral *> gateName2Gate;

    // introduce a Node for each variable
    for (auto var : res.freeVariables)
        gateName2Gate[var] = Variable::getVariableGate(var);
    for (auto e : res.existentialVariables)
        for (auto var : e.second)
            gateName2Gate[var] = Variable::getVariableGate(var);
    for (auto u : res.universalVariables)
        for (auto var : u.second)
            gateName2Gate[var] = Variable::getVariableGate(var);

    // read lines until end of file
    while (std::getline(input, line))
    {
        // printf("Line: %s\n", line.c_str());
        if (line[0] == '#') // ignore comments
            continue;

        // split at '=' sign
        std::size_t found = line.find("=");
        if (found == std::string::npos)
        {
            printf("Error: Line doesn't contain gate defintion\n");
            EXIT_UNWANTED_STATE
        }
        int gateName = stoi(trim(line.substr(0, found)));
        std::string gateDef = trim(line.substr(found + 1)); // ignore opening bracket

        std::vector<int> inputs;
        if (parseLinePrefix(gateDef, "or", inputs))
        {
            vector<signedGate_t> inputGates;
            for (auto input : inputs)
            {
                int sign = input < 0 ? -1 : 1;
                inputGates.push_back({sign, gateName2Gate[abs(input)]});
            }
            OrGate *gate = OrGate::getOrGate(inputGates);
            gateName2Gate[gateName] = gate;
        }
        else if (parseLinePrefix(gateDef, "and", inputs))
        {
            vector<signedGate_t> inputGates;
            for (auto input : inputs)
            {
                int sign = input < 0 ? -1 : 1;
                inputGates.push_back({sign, gateName2Gate[abs(input)]});
            }
            AndGate *gate = AndGate::getAndGate(inputGates);
            gateName2Gate[gateName] = gate;
        }
        else
        {
            printf("Error: Unknown gate type: %s\n", line.c_str());
            EXIT_UNWANTED_STATE
        }
    }

    res.outputGate = gateName2Gate[outputGateName];
    return res;
}

// check whether the input is a valid 2QBF instance
void checkInput(qcir_t &instance)
{
    if (instance.existentialVariables.size() > 1)
    {
        printf("Error: Only one existential level is allowed\n");
        EXIT_UNWANTED_STATE
    }

    if (instance.universalVariables.size() > 1)
    {
        printf("Error: Only one universal part is allowed which must come after the existential one\n");
        EXIT_UNWANTED_STATE
    }

    // if (instance.universalVariables.size() == 0)
    // {
    //     printf("Error: No univeral part present\n");
    //     EXIT_UNWANTED_STATE
    // }

    if (instance.existentialVariables.size() == 1 && instance.universalVariables.size() == 1)
    {
        if (instance.existentialVariables[0].first > instance.universalVariables[0].first)
        {
            printf("Error: Universal part must come after the existential one\n");
            EXIT_UNWANTED_STATE
        }
    }
}