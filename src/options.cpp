#include "options.hpp"

po::options_description main_opts("Main options");
po::options_description output_opts("Output options");
po::options_description minimality_opts("Minimality checker options");
po::options_description propagator_opts("Propagator options");
po::options_description perf_opts("Performance options");
po::options_description other_opts("Other options");

po::options_description all_opts("All options");

void initOptions(SolverConfig &config, struct minimality_config_t &minimalityConfig, propagators_config_t &propagatorsConfig)
{
    // clang-format off
    main_opts.add_options()
      ("all-graphs,a", po::bool_switch(&config.allModels), "Generate all graphs until unsatisfiable (otherwise exit after first solution)")
      ("vertices,v", po::value<int>(&config.vertices), "Generate graphs with this number of vertices")
      ("directed,d", po::bool_switch(&config.directed), "Generate directed graphs")
      ("no-SMS", po::bool_switch(&minimalityConfig.turnoffSMS), "Turn off SMS, i.e., no minimality check")
      ("dimacs", po::value<string>() ,"File containing the encoding for the wished for graphs, in DIMACS format")
      ;

    output_opts.add_options()
      ("help", "Produce this help message")
      ("hide-graphs", po::bool_switch(&config.hideGraphs), "Count solutions without outputing them (they still have to be enumerated)")
      ;

    minimality_opts.add_options()
      ("frequency", po::value<int>(&minimalityConfig.frequency)->default_value(DEFAULT_MINIMALITY_FREQUENCY), "The inverse frequency with which to call the minimality check, i.e., call every f-th time the propagator is called")
      ("cutoff", po::value<int>(&minimalityConfig.cutoff)->default_value(DEFAULT_MINIMALITY_CUTOFF), "The cutoff for the minimality check (0 for deactivating the cutoff).")
      ("initial-partition", po::value<vector<int>>()->multitoken(), "Set the initial partition for the minimality check, given as a sequence of partition sizes (integers). The minimality check will look for permutations that preserve partition membership")
      ;
      
    // try to sort them by general relevance
    propagator_opts.add_options()
      ("connected,c", po::bool_switch(&propagatorsConfig.generate_connected), "Generate only connected graphs")
      ("planar,p", po::bool_switch(&propagatorsConfig.planar), "Generate planar graphs")
      ("planarity-frequency", po::value<int>(&propagatorsConfig.planarityFrequency)->default_value(DEFAULT_FREQUENCY), "The frequency with which to call the planarity check")
      ("forbidden-subgraph-file", po::value<string>(&propagatorsConfig.forbiddenSubgraphFile), "File containing a list of forbidden subgraphs")
      ("forbidden-induced-subgraph-file", po::value<string>(&propagatorsConfig.forbiddenInducedSubgraphFile), "File containing a list of forbidden induced subgraphs")
      ("frequency-forbidden-subgraphs", po::value<int>(&propagatorsConfig.frequencyForbiddenSubgraphs), "The frequency with which to call the forbidden (induced) subgraph check")
      ("min-chromatic-number,chi", po::value<int>(&propagatorsConfig.minChromaticNumber)->implicit_value(0), "Search for graphs with chromatic number at least this")
      ("coloring-algo", po::value<int>(&propagatorsConfig.coloringAlgo), "Specify the graph coloring algorithm (0 means simple recursive; 1 means simple DPLL-based; 2 means SAT-based)")
      ("non-010-colorable,non010", po::bool_switch(&propagatorsConfig.non010colorable), "Search for non-010-colorable graphs (relates to Kochen-Specker graphs)")
      ("triangle-vars", po::value<int>(&propagatorsConfig.triangleVars)->implicit_value(0), "For problems which use triangle variables (variable to denote whether a given set of 3 vertices forms a triangle), specify where the triangle variables start")
      ;
    

    other_opts.add_options()
      ("simplify", po::value<std::string>(&config.simplifiedFormulaFile), "Simplify the CNF formula and write it to the given file. (Directly after prerun)")
      ("learned-clauses", po::value<std::string>(&config.learnedClausesFile),"Write the learned clauses to the given file after prerun")
      ("max-learned-clause-size", po::value<int>(&config.maxPrintedLearnedClauseSize)->default_value(5), "The maximal size of the learned clauses that should be printed (also for `simplify`)")
      ("prerun", po::value<int>(&config.prerunTime), "Run the solver for the given amount of seconds before creating cubes/ simplified formula/ output learned clauses/ activating lookahead.")
      ("assignment-cutoff", po::value<int>(&config.assignmentCutoff), "Create cubes when the number of assigned variables exceeds this (using lookahead heuristic for branching)") // TODO distinguish between different cubing types
      ("simple-assignment-cutoff", po::value<int>(&config.simpleAssignmentCutoff), "Create cubes soley based on edge variables. Cubes are blocked by clauses until unsat is reached")
      ("create-game", po::bool_switch(&config.createGame), "Print the trail on a conflict or if all clauses are satisfied using lookeahead")
      ("create-randomized-game", po::value<float>(&config.createGameProp), "Randomize the truth value of the branching variable. The value is the probability of flipping the truth value.")
      ("create-game-recursion-lvl", po::value<int>(&config.createGameRecLvl), "The depth of the lookahead on the truth values for creating games")
      ("cube-file", po::value<std::string>(&config.cubeFile), "Filename containing the cubes")
      ("cube-file-test", po::value<std::string>(&config.cubeFileTest), "Test the given cubes, i.e., whether negation of the cubes is unsat")
      ("cube-line", po::value<int>(&config.cubeLine), "The cube from the given line is selected (starting with line 1)")
      ("cubes-range", po::value<vector<int>>(&config.cubesRange)->multitoken(), "Range of cubes to be processed")
      ("cube-timeout", po::value<int>(&config.cubeTimeout), "Maximum time spend for trying to solve for each cube in seconds")
      ("timeout", po::value<int>(&config.timeout), "Timeout in seconds of the main solving phase")
      ("cube-only-decisions",po::bool_switch(&config.cubesOnlyDecisions), "The output of the cubes only contains the decision literals and not implied literals")
      ("lookahead-only-edge-vars", po::bool_switch(&config.lookaheadOnlyEdgeVars), "Use lookahead only on edge variables (either for cubing of creating games)")
      ("cadical-config", po::value<vector<std::string>>(&config.cadicalConfig)->multitoken(), "Space-separated list of command-line options for cadical (without the -- prefixes)")
      ("qcir", po::value<string>(&propagatorsConfig.qcirFile) ,"File containing the encoding in QCIR format. (Currently only supports 2QBF starting with existential quantifier)")
      ("polarity-hasing", po::bool_switch(&propagatorsConfig.polarityHashing), "Use polarity hashing for the 2QBF solver")
      ;
    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------

    perf_opts.add_options(); // TODO
    
    all_opts.add(main_opts).add(output_opts).add(minimality_opts).add(propagator_opts).add(perf_opts).add(other_opts);
    // clang-format on
}

#include "graphPropagators/connectedChecker.hpp"
#include "graphPropagators/planarity.hpp"
#include "graphPropagators/coloringCheck.hpp"
#include "other/forbiddenSubgraph.hpp"
#include "qbf/universal2.hpp"

void addPropagators(GraphSolver *solver, const propagators_config_t &propagatorsConfig)
{
    vector<vector<int>> edgeVariables = solver->getEdgeVariables();

    if (propagatorsConfig.generate_connected)
        solver->addPartiallyDefinedGraphChecker(new ConnectedChecker());

    if (propagatorsConfig.planar)
    {
        if (solver->config.directed)
            solver->addPartiallyDefinedGraphChecker(new DirectedPlanarityChecker(propagatorsConfig.planarityFrequency));
        else
            solver->addPartiallyDefinedGraphChecker(new PlanarityChecker(propagatorsConfig.planarityFrequency));
    }

    if (propagatorsConfig.minChromaticNumber > 0)
        solver->addComplexFullyDefinedGraphChecker(new MinChromaticNumberChecker(propagatorsConfig.minChromaticNumber, propagatorsConfig.coloringAlgo, edgeVariables, false));

    if (propagatorsConfig.non010colorable)
    {
        int n = solver->config.vertices;
        vector<vector<vector<int>>> triangleVars(n, vector<vector<int>>(n, vector<int>(n, 0)));

        if (!propagatorsConfig.triangleVars)
        {
            printf("Starting variable for triangle vars must be specified when searching for non-010-colorable graphs\n");
            EXIT_UNWANTED_STATE
        }

        int c = propagatorsConfig.triangleVars;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                for (int k = j + 1; k < n; k++)
                    triangleVars[i][j][k] = c++;

        solver->addComplexFullyDefinedGraphChecker(new Non010colorableChecker(triangleVars, edgeVariables));
    }

    if (propagatorsConfig.forbiddenSubgraphFile != "" || propagatorsConfig.forbiddenInducedSubgraphFile != "")
    {
#ifdef GLASGOW
        std::ifstream forbiddenSubgraphFile;
        std::ifstream forbiddenInducedSubgraphFile;
        std::ifstream forbiddenPartialSubgraphFile;

        if (propagatorsConfig.forbiddenSubgraphFile != "")
            forbiddenSubgraphFile.open(propagatorsConfig.forbiddenSubgraphFile);
        if (propagatorsConfig.forbiddenInducedSubgraphFile != "")
            forbiddenInducedSubgraphFile.open(propagatorsConfig.forbiddenInducedSubgraphFile);
        solver->addPartiallyDefinedGraphChecker(new ForbiddenSubgraphCheckerGlasgow(propagatorsConfig.frequencyForbiddenSubgraphs, forbiddenSubgraphFile, forbiddenInducedSubgraphFile, forbiddenPartialSubgraphFile));
#else
        LOG(LOG_LEVEL_INFO, "Forbidden subgraph check is not available. Recompile with the appropriate flag.");
        throw std::runtime_error("Forbidden subgraph check is not available. Recompile with the appropriate flag.");
#endif
    }

    if (propagatorsConfig.qcirFile != "")
    {
        std::ifstream qcirFile(propagatorsConfig.qcirFile);
        if (!qcirFile.is_open())
        {
            throw std::runtime_error("Could not open QCIR file");
        }
        QCIRchecker *qcirChecker = new QCIRchecker(qcirFile, propagatorsConfig.polarityHashing);
        solver->addComplexFullyDefinedGraphChecker(qcirChecker);
        // set highest variable !!!
        solver->numVars = std::max(solver->numVars, qcirChecker->highestVariableInInstance);
        printf("Highest variable in instance: %d\n", qcirChecker->highestVariableInInstance);
        // add dummy clause so that the solver knows the highest variable
        solver->add(qcirChecker->highestVariableInInstance);
        solver->add(-qcirChecker->highestVariableInInstance);
        solver->add(0);
    }
}

GraphSolver *createSolver(SolverConfig &config, struct minimality_config_t &minimalitConfig, propagators_config_t &propagatorsConfig, string dimacsFile)
{
    GraphSolver *solver = new GraphSolver(config, minimalitConfig);
    solver->set("quiet", 1);

    if (!dimacsFile.empty())
    {
        int maxVar;
        solver->read_dimacs(dimacsFile.c_str(), maxVar);
        solver->numVars = std::max(solver->numVars, maxVar);
    }

    if (!config.cubeFileTest.empty())
    {
        // add negation of cubes as clauses
        std::ifstream cubeFile(config.cubeFileTest);
        if (!cubeFile.is_open())
        {
            throw std::runtime_error("Could not open cube file");
        }

        std::string line;
        while (std::getline(cubeFile, line))
        {
            // skip lines not starting with "a "
            if (line.substr(0, 2) != "a ")
                continue;
            std::istringstream iss(line);
            std::string lit;
            clause_t clause;
            while (iss >> lit)
            {
                if (lit == "a")
                    continue;
                int litInt = std::stoi(lit);
                solver->add(-litInt); // add negation
            }
        }
    }

    addPropagators(solver, propagatorsConfig);

    // --------- check some invalid or problematic combinations ------------
    if (config.allModels && minimalitConfig.cutoff != 0)
        LOG(LOG_LEVEL_INFO, "The output might contain isomorphic copies of the same graph. Set 'cutoff' to 0 to avoid this potentially resulting in a performance loss.");
    return solver;
}

GraphSolver *parseOptions(const int argc, const char **argv, po::variables_map &vm)
{
    SolverConfig config;
    struct minimality_config_t minimalitConfig;
    propagators_config_t propagatorsConfig;

    initOptions(config, minimalitConfig, propagatorsConfig);

    try
    {
        po::store(po::parse_command_line(argc, argv, all_opts), vm);
        po::notify(vm);
    }
    catch (const po::unknown_option &e)
    {
        std::cout << "Error: Unknown option '" << e.get_option_name() << "'" << std::endl;
        std::cout << "Note that when using pysms the unknown arguments are forwarded to smsg/smsd." << std::endl;
        std::cout << all_opts << std::endl;
        exit(1);
    }

    if (vm.count("help"))
    {
        std::cout << all_opts << std::endl;
        exit(0);
    }

    // create the solver

    string dimacsFile;
    if (vm.count("dimacs"))
        dimacsFile = vm["dimacs"].as<string>();

    if (vm.count("initial-partition"))
    {
        vector<int> partitionSizes = vm["initial-partition"].as<vector<int>>();
        minimalitConfig.initialPartition = vector<bool>(config.vertices, false);
        int curPos = 0;
        for (const int value : partitionSizes)
        {
            minimalitConfig.initialPartition[curPos] = true;
            curPos += value;
        }
        assert(config.vertices == curPos);
    }

    GraphSolver *solver = createSolver(config, minimalitConfig, propagatorsConfig, dimacsFile);
    return solver;
}
