""" 
    A script for creating graph encodings in QCIR format.
"""

from itertools import *
from sys import *
import argparse
import os


class IDPool:
    """A class for returning the next free id starting with 1"""

    def __init__(self, start_from=1) -> None:
        self.nextId = start_from

    def id(self):
        # returns the next free variable
        x = self.nextId
        self.nextId += 1
        return x


class Gate:
    def __init__(self, id, inputs=None) -> None:
        self._inputs = inputs
        self._id = id
        self._gateName = ""

    def __getInputsAsString(self):
        return ",".join(map(lambda x: str(x) if type(x) == int else str(x.getSignedId()), self._inputs))

    def getId(self):
        return self._id

    def getSignedId(self):
        return self.getId()

    def __neg__(self):
        return NegatedGate(self)

    def gateEncoding(self):
        return f"{self._id} = {self.gateName}({self.__getInputsAsString()})"

    def setInput(self, newInput):
        self._inputs = newInput

    def appendToInput(self, lit):
        self._inputs.append(lit)

    def getInput(self):
        return self._inputs


class NegatedGate(Gate):
    """Takes gate but negation"""

    def __init__(self, gate: Gate) -> None:
        self.__gate = gate

    def getId(self):
        return self.__gate.getId()

    def getSignedId(self):
        return -self.__gate.getSignedId()

    def __neg__(self):
        return self.__gate  # double negation is the oridignal gate

    def gateEncoding(self):
        return self.__gate.gateEncoding()

    def setInput(self, newInput):
        self.__gate.setInput(newInput)

    def appendToInput(self, lit):
        self.__gate.appendToInput(lit)

    def getInput(self):
        return self.__gate.getInput()


class AndGate(Gate):
    def __init__(self, id, inputs=None) -> None:
        super().__init__(id, inputs)
        self.gateName = "and"


class OrGate(Gate):
    def __init__(self, id, inputs=None) -> None:
        super().__init__(id, inputs)
        self.gateName = "or"


def seqCounter(variables, countUpto, vPool, outputGate: Gate, atMost=None, atLeast=None):
    x = vPool.id()
    FalseGate = AndGate(vPool.id(), [x, -x])  # for the sake of simplicity

    n = len(variables)
    counterVariables = [[OrGate(vPool.id(), []) for _ in range(countUpto)] for _ in range(n)]  # Create new gates
    # print("c\t" + str(counterVariables))
    # first element
    counterVariables[0][0] = variables[0]
    for i in range(1, countUpto):
        counterVariables[0][i] = FalseGate  # at most one at the beginning

    # adapt counter for each step
    for i in range(n - 1):
        counterVariables[i + 1][0].appendToInput(variables[i + 1])  # if there is an element than there is at least on element
        for j in range(countUpto):
            counterVariables[i + 1][j].appendToInput(counterVariables[i][j])  # at least as many

            # clauses.append([counterVariables[i][j], variables[i + 1], -counterVariables[i + 1][j]])  # the same if element is not present

            if j < countUpto - 1:
                counterVariables[i + 1][j + 1].appendToInput(AndGate(vPool.id(), [counterVariables[i][j], variables[i + 1]]))  # one more element
                # clauses.append([counterVariables[i][j], -counterVariables[i + 1][j + 1]])  # at most one more

    if atMost:
        for i in range(n - 1):
            g = OrGate(vPool.id(), [-counterVariables[i][atMost - 1], -variables[i + 1]])
            outputGate.appendToInput(g)  # if maximum reached, no more true variables
    if atLeast:
        outputGate.appendToInput(OrGate(vPool.id(), [counterVariables[n - 1][atLeast - 1]]))

    return [counterVariables[n - 1][j] for j in range(countUpto)]


def getDefaultParser():
    parser = argparse.ArgumentParser()  # default parser for some properties

    main_args = parser.add_argument_group(title="Main arguments", description="The number of vertices is mandatory, everything else is optional")
    solve_args = parser.add_argument_group(title="Solver options")
    constraint_args = parser.add_argument_group(
        title="Graph constraints", description="A set of pre-defined constraints for common applications, including applications from SMS papers"
    )

    main_args.add_argument("--vertices", "-v", type=int, required=True, help="number of vertices")
    main_args.add_argument("--directed", "->", action="store_true", help="search for directed graphs")
    main_args.add_argument("--print-dimacs", action="store_true", help="Print output in DIMACS format assuming no quantifiers")

    # solve options
    solve_args.add_argument(
        "--all-graphs", "-a", action="store_true", help="generate all graphs (without this, the solver will exit after the first solution)"
    )
    solve_args.add_argument(
        "--hide-graphs",
        "-hg",
        action="store_true",
        help="do not display graphs (meant as a counting functionality, though the graphs still need to be enumerated)",
    )
    solve_args.add_argument(
        "--args-SMS", type=str, default="", help="command line to be appended to the call to smsg/smsd (see src/main.cpp or README.md)"
    )
    solve_args.add_argument(
        "--print-cnf",
        type=str,
        help="don't run SMS, only generate a CNF encoding of the gate omitting the quantifiers and print it to the given file",
    )
    solve_args.add_argument("--print-qcir", type=str, help="don't run SMS, only generate a QCIR encoding of the gate and print it to the given file")
    solve_args.add_argument("--qcir-file", type=str, help="file to store the temporary qcir encoding later used by SMS")
    solve_args.add_argument("--graph6-format", action="store_true", help="output graphs in graph6 format (Warning: Relatively slow)")

    constraint_args.add_argument("--negate-output-gate", action="store_true", help="negate the output gate (but doesn't swap quantifiers)")

    # number of edges
    constraint_args.add_argument("--num-edges-upp", type=int, help="upper bound on the maximum number of edges")
    constraint_args.add_argument("--num-edges-low", type=int, help="lower bound on the minimum number of edges")

    # degree restrictions
    constraint_args.add_argument("--Delta-upp", type=int, help="upper bound on the maximum degree")
    constraint_args.add_argument("--delta-low", type=int, help="lower bound on the minimum degree")
    constraint_args.add_argument("--cubic", action="store_true", help="ensure that the graph is cubic, i.e., 3-regular")

    # cycles
    constraint_args.add_argument("--Ck-free", type=int, help="forbid the k-cycle C_k as (non-induced) subgraph")

    # treewidth
    constraint_args.add_argument(
        "--tree-width-upp-version1",
        type=int,
        help="ensure that the treewidth is at most a certain value using a compact encoding but using more variables",
    )
    constraint_args.add_argument("--tree-width-upp-version2", "--tw", type=int, help="ensure that the treewidth is at most a certain value")
    constraint_args.add_argument(
        "--tree-width-low-version1",
        type=int,
        help="ensure that the treewidth is at least a certain value using a compact encoding but using more variables",
    )
    constraint_args.add_argument("--tree-width-low-version2", type=int, help="ensure that the treewidth is at least a certain value")

    constraint_args.add_argument(
        "--tree-width-vertex-critical", "--twvc", type=int, help="ensure that deleting any vertex results in treewidth the given value - 1"
    )
    constraint_args.add_argument("--tree-width-critical", "--twc", type=int, help="ensure that every minor has treewidth given value - 1")

    # coloring
    constraint_args.add_argument("--chi-low", type=int, help="Minimum chromatic number")
    constraint_args.add_argument("--chi-upp", type=int, help="Maximum chromatic number")

    constraint_args.add_argument("--non-010", action="store_true", help="ensure that non 010 colorable")
    constraint_args.add_argument("--kochen-specker", action="store_true", help="ensure that candidate for kochen specker graph")

    constraint_args.add_argument("--non-3-edge-colorable", action="store_true", help="ensure that non 3 edge colorable")

    # side constraints
    constraint_args.add_argument("--no-subsuming-neighborhoods", action="store_true", help="ensure that no subsuming neighborhoods are present")
    constraint_args.add_argument("--mtf", action="store_true", help="ensure that maximal triangle free")
    constraint_args.add_argument("--connected-static", action="store_true", help="Static encoding for being connected")
    constraint_args.add_argument("--two-connected", action="store_true", help="deleting any vertex doesn't make the graph disconnected")

    constraint_args.add_argument(
        "--folkmann",
        type=int,
        help="ensures that each 2 coloring has a monochromatic edge and that there is no clique of size given value (some additional constraints are added)",
    )

    constraint_args.add_argument("--maximal-4-clique-free", action="store_true", help="no 4-clique and inserting any edge results in a 4-clique")

    constraint_args.add_argument("--max-independent-set", type=int, help="maximum independent set of given size")

    constraint_args.add_argument("--min-domination-number", type=int, help="lower bound on the domation number")
    constraint_args.add_argument(
        "--three-decomposition-conjecture", action="store_true", help="search for counter example to 3-decomposition conjecture"
    )
    constraint_args.add_argument("--bipartite", action="store_true", help="ensure that graph is bipartite")
    constraint_args.add_argument(
        "--domination-conjecture", action="store_true", help="ensure that resulting graph is cubic and has domination number > ceil(n/3) "
    )
    constraint_args.add_argument(
        "--min-girth-compact", type=int, help="compact encoding to ensure that the girth is at least k, i.e., no cycle with length < k"
    )

    constraint_args.add_argument("--canonical-qbf", action="store_true", help="use universal variables to ennsure that the graph is lex-minimal")
    constraint_args.add_argument("--canonical-qbf-colex", action="store_true", help="same as canonical-qbf but using colex ordering")
    return parser


DEFAULT_COUNTER = "sequential"


def transversal_rec(curGate, alreadyExpanded):
    """Function transversing all gates and variables starting with a specific gate"""
    # print(type(curGate))
    if not isinstance(curGate, Gate):
        yield abs(curGate)
    else:
        for g in curGate.getInput():
            gate_id = g.getId() if isinstance(g, Gate) else abs(g)
            if not gate_id in alreadyExpanded:  # if already expanded then done
                alreadyExpanded.add(gate_id)
                yield from transversal_rec(g, alreadyExpanded)
        if isinstance(curGate, NegatedGate):
            yield -curGate
        else:
            yield curGate


class GraphEncodingBuilder(IDPool):
    """A class for building an encoding in the context of SMS.
    The IDPool gives the next free variable"""

    def __init__(self, n, directed=False):
        super().__init__()
        self.directed = directed
        self.V = list(range(n))
        self.n = n
        self.paramsSMS = {"vertices": self.n, "print-stats": True, "frequency": 30, "cutoff": 20000}
        self.varEdgeTable = [[None for _ in self.V] for _ in self.V]
        self.varEdgeDirectedTable = [[None for _ in self.V] for _ in self.V]
        if directed:
            for v, u in permutations(self.V, 2):
                self.varEdgeDirectedTable[v][u] = self.id()
        for v, u in combinations(self.V, 2):
            self.varEdgeTable[v][u] = self.varEdgeTable[u][v] = self.id()

        self.outputGate = AndGate(self.id(), [])  # and gate collecting all the constraints
        self.quantifiers = [
            ("exists", []),
            ("forall", []),
        ]  # list of tuples containing the quantification type and the varialbes of the level for example [("exits", [10,42]), ]. Can be expanded by further quantifier levels
        self.DEBUG = True

    def var_edge(self, u, v) -> int:
        """Get propositional variable associated with the undirected edge u,v"""
        return self.varEdgeTable[u][v]

    def var_edge_dir(self, u, v) -> int:
        """Get propositional variable associated with the directed edge (u,v)"""
        return self.varEdgeDirectedTable[u][v]

    # ------------------some utilities--------------------------

    def sorted(self):
        yield from transversal_rec(self.outputGate, set())

    def transversal_rec(curGate, alreadyExpanded):
        """Function transversing all gates and variables (and returning them) starting with a specific gate"""
        # print(type(curGate))
        if not isinstance(curGate, Gate):
            yield abs(curGate)
        else:
            for g in curGate.getInput():
                gate_id = g.getId() if isinstance(g, Gate) else abs(g)
                if not gate_id in alreadyExpanded:  # if already expanded then done
                    alreadyExpanded.add(gate_id)
                    yield from transversal_rec(g, alreadyExpanded)
            if isinstance(curGate, NegatedGate):
                yield -curGate
            else:
                yield curGate

    def getVariables(gate):
        gatesAndVariables = transversal_rec(gate, set())
        variables = []
        for x in gatesAndVariables:
            if isinstance(x, int):
                variables.append(x)
        return variables

    def getFreeVariables(self):
        """Returns all edge and directed edge variables"""
        return (set(chain.from_iterable(self.varEdgeTable)) | set(chain.from_iterable(self.varEdgeDirectedTable))) - {None}

    def addExistentialGate(self, g):
        """
        Add a gate to the existential part of the formula.
        """
        self.outputGate.appendToInput(g)
        # add variables to forall block
        assert self.quantifiers[0][0] == "exists"
        edge_vars = self.getFreeVariables()
        self.quantifiers[0][1].extend([v for v in GraphEncodingBuilder.getVariables(g) if v not in edge_vars])

    def addUniversalGate(self, g):
        """
        Add a gate to the universal part of the formula.
        All variables except edge variables are added to the universal block
        """
        self.outputGate.appendToInput(g)
        # add variables to forall block
        assert self.quantifiers[1][0] == "forall"
        edge_vars = self.getFreeVariables()
        self.quantifiers[1][1].extend([v for v in GraphEncodingBuilder.getVariables(g) if v not in edge_vars])

    def print_qcir(self, outputFile=stdout):
        """Prints the output gate in QCIR format.
        The quantifiers are specified by self.quantifiers.
        Edge variables are free.
        """
        print("#QCIR-G14", file=outputFile)
        inputs = []
        gates = []
        for g in self.sorted():
            if isinstance(g, Gate):
                # print("Gate:", g.getId())
                gates.append(g)
            else:
                # print("Var:", g)
                inputs.append(g)

        edge_vars = self.getFreeVariables()
        free_inputs = sorted(edge_vars)

        print(f"free({','.join(map(str, free_inputs))})", file=outputFile)
        #
        for t, variables in self.quantifiers:
            if len(variables) == 0:
                continue
            print(f"{t}({','.join(map(str, variables))})", file=outputFile)
        print(f"output({self.outputGate.getId()})", file=outputFile)

        for g in gates:
            print(g.gateEncoding(), file=outputFile)

    def print_dimacs(self, output=stdout):
        clauses = [[self.outputGate.getId()]]
        # output must be true
        # assume that no quantifiers and write as cnf using teitsin transformation
        # TODO polarities and further optimizations
        for g in self.sorted():
            if not isinstance(g, Gate):
                continue

            if isinstance(g, OrGate):
                clauses.append([-g.getId()] + [h.getSignedId() if isinstance(h, Gate) else h for h in g.getInput()])
                for h in g.getInput():
                    clauses.append([g.getId(), -h.getSignedId() if isinstance(h, Gate) else -h])

            elif isinstance(g, AndGate):
                clauses.append([g.getId()] + [-h.getSignedId() if isinstance(h, Gate) else -h for h in g.getInput()])
                for h in g.getInput():
                    clauses.append([-g.getId(), h.getSignedId() if isinstance(h, Gate) else h])
            else:
                print("Error: current gate type not supported", type(g), file=stderr)

        maxVariable = max(abs(value) for row in clauses for value in row)
        print(f"p cnf {maxVariable} {len(clauses)}", file=output)
        for c in clauses:
            print(" ".join(str(x) for x in c), 0, file=output)

    # ----------- list of constraints without quantification, must be added manually

    def colorable(self, chi):
        # ensure that graph is colorable with chi colors
        outputGate = AndGate(self.id(), [])
        coloring = [[self.id() for c in range(chi)] for v in self.V]

        # at least one color
        for v in self.V:
            outputGate.appendToInput(OrGate(self.id(), coloring[v]))

        for u, v in combinations(self.V, 2):
            for c in range(chi):
                outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v), -coloring[u][c], -coloring[v][c]]))
        return outputGate

    def colorable010(self, triangleVars=None):
        # ensure that graph is colorable with chi colors
        outputGate = AndGate(self.id(), [])
        coloring = [self.id() for v in self.V]

        if triangleVars == None:
            print("Error: triangleVars must be given")
            exit(1)
            # triangleVars = {(u, v, w): AndGate(self.id(), [self.var_edge(u, v), self.var_edge(v, w), self.var_edge(u, w)]) for u, v, w in combinations(self.V, 3)}

        for u, v in combinations(self.V, 2):
            outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v), -coloring[u], -coloring[v]]))

        for u, v, w in combinations(self.V, 3):
            triangleVars[(u, v, w)]
            outputGate.appendToInput(OrGate(self.id(), [-triangleVars[(u, v, w)], coloring[u], coloring[v], coloring[w]]))
        return outputGate

    def dominatingSet(self, k):
        """Ensure that graph has a dominating set of size at most k"""
        outputGate = AndGate(self.id(), [])
        selectedVertices = [self.id() for v in self.V]
        # at most k selected vertices
        seqCounter(selectedVertices, k, self, outputGate, atMost=k)

        for v in self.V:
            # either selcted or neighbor selected
            g = OrGate(self.id(), [+selectedVertices[v]] + [AndGate(self.id(), [self.var_edge(v, u), selectedVertices[u]]) for u in self.V if u != v])
            outputGate.appendToInput(g)
        return outputGate

    def minGirthCompact(self, minGirth):
        """Compact encoding to ensure that the girth is at least k, i.e., no cycle with length < k
        :param minGirth: girth lower bound
        """
        assert minGirth >= 4

        # check distance of i,j without edge i,j.
        for i, j in combinations(self.V, 2):

            reached = [self.var_edge(i, k) if k not in [i, j] and k > i else None for k in self.V]
            # the outcommented lines are the CNF version of the encoding
            for _ in range(minGirth - 4):  # if girth is 4 than no triangles so not in the loop
                # reachedNew = [self.id() for _ in self.V]
                reachedNew = [OrGate(self.id(), []) for _ in self.V]
                for k in self.V:
                    if k in [i, j] or k < i:
                        continue
                    # self.append([-reached[k], +reachedNew[k]])  # already reached
                    reachedNew[k].appendToInput(reached[k])

                    # check if reached over l
                    for l in self.V:
                        if l in [i, j, k] or l < i:
                            continue
                        # self.append([-self.var_edge(k, l), -reached[l], +reachedNew[k]])  # l reached in previous step and edge implies reached
                        reachedNew[k].appendToInput(AndGate(self.id(), [self.var_edge(k, l), reached[l]]))
                reached = reachedNew

            for k in self.V:
                if k in [i, j] or k < i:
                    continue
                # not reached, not adjacent to j, or edge not present
                # self.append([-self.var_edge(i, j), -self.var_edge(j, k), -reached[k]])
                self.outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(i, j), -self.var_edge(j, k), -reached[k]]))

    def lexSmallerMatrix(self, useColexOrdering):
        """Ensure that there is a permutation which makes the adjacency matrix lexicographically smaller, i.e., there is an indicator pair"""
        assert not self.directed
        outputGate = AndGate(self.id(), [])

        p = [[self.id() for _ in self.V] for _ in self.V]  # permutation matrix
        # ensure that permutation i.e., total and injective
        for i in self.V:
            outputGate.appendToInput(OrGate(self.id(), p[i]))  # at least mapped to one position
        for i, j in combinations(self.V, 2):
            for k in self.V:
                outputGate.appendToInput(OrGate(self.id(), [-p[i][k], -p[j][k]]))  # injective

        # get all vertex pairs mapped to the same vertex
        # precEqual = [[OrGate(self.id(), []) if i < j else None for j in self.V] for i in self.V]  # whether edge i,j is mapped to edge i,j
        # for i, j in combinations(self.V, 2):
        #     precEqual[i][j].appendToInput(AndGate(self.id(), [p[i][i], p[j][j]]))
        #     precEqual[i][j].appendToInput(AndGate(self.id(), [p[i][j], p[j][i]]))

        # get permuted matrix
        matrixPermuted = [[OrGate(self.id(), []) if i < j else None for _ in self.V] for _ in self.V]  # permuted matrix
        for i, j in combinations(self.V, 2):
            # edge between i,j iff there is u,v such that edge u,v is mapped to i,j
            for u, v in permutations(self.V, 2):
                matrixPermuted[i][j].appendToInput(AndGate(self.id(), [self.var_edge(u, v), p[i][u], p[j][v]]))

        # there must be at least one indicator pair
        indicatorPairPresent = OrGate(self.id(), [])

        # if all previous are greater equal then originial matrix either 0 or permuted one 1
        allPrevGreaterEqual = AndGate(self.id(), [])  # empty and get is true
        vertexOrdering = list(combinations(self.V, 2))
        if useColexOrdering:
            vertexOrdering = []  # colex ordering
            for j in self.V:
                for i in self.V:
                    if i < j:
                        vertexOrdering.append((i, j))
        for i, j in vertexOrdering:
            indicatorPairPresent.appendToInput(AndGate(self.id(), [allPrevGreaterEqual, self.var_edge(i, j), -matrixPermuted[i][j]]))
            greaterEqual = OrGate(
                self.id(), [self.var_edge(i, j), -matrixPermuted[i][j]]
            )  # , precEqual[i][j]])  # check whether i,j is greater equal
            allPrevGreaterEqual = AndGate(self.id(), [allPrevGreaterEqual, greaterEqual])  # all previous equal and new either equal or an edge

        outputGate.appendToInput(indicatorPairPresent)
        # outputGate.appendToInput(OrGate(self.id(), []))  # False gate for testing
        return outputGate

    # check whether the graph given by the edge variabes/gates is connected
    def connectedStatic(self, n, var_edge=None):
        if var_edge == None:
            var_edge = self.var_edge

        outputGate = AndGate(self.id(), [])
        startingVertex = 0
        R = list(range(1, n))  # all vertices accept the first

        reached = {v: var_edge(startingVertex, v) for v in R}
        for k in reached:

            assert reached[k] != None
        reachedNew = {v: OrGate(self.id(), []) for v in R}
        for _ in R:  # number of rounds
            for v in R:
                reachedNew[v].appendToInput(reached[v])
            for v, w in permutations(R, 2):
                reachedNew[v].appendToInput(AndGate(self.id(), [reached[w], var_edge(v, w)]))
            reached = reachedNew
            reachedNew = {v: OrGate(self.id(), []) for v in R}

        for v in R:
            outputGate.appendToInput(reached[v])
        return outputGate

    def bipartite(self):
        """Ensure that graph is bipartite"""
        outputGate = AndGate(self.id(), [])
        color1 = [self.id() for v in self.V]  # either in one or the other partition
        for u, v in combinations(self.V, 2):
            outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v), -color1[u], -color1[v]]))
            outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v), +color1[u], +color1[v]]))
        outputGate.appendToInput(OrGate(self.id(), [color1[0]]))  # fix first vertex
        return outputGate

    def specialSpanningTree(self):
        """Ensure that there is a spanning tree such that the remaining graph can be split into a mathing and a 2-regular graph (see 3-Decomposition Conjecture)"""
        # No to remaining vertices where one of them has degree 1
        outputGate = AndGate(self.id(), [])
        spanningTreeEdges = [[self.id() if u < v else None for v in self.V] for u in self.V]

        def var_span_edge(u, v):
            return spanningTreeEdges[min(u, v)][max(u, v)]

        for u, v in combinations(self.V, 2):
            outputGate.appendToInput(OrGate(self.id(), [-var_span_edge(u, v), self.var_edge(u, v)]))

        # ensure that it is a spanning tree
        numEdges = self.n - 1
        seqCounter([var_span_edge(u, v) for u, v in combinations(self.V, 2)], numEdges, self, outputGate, atLeast=numEdges, atMost=numEdges)

        # ensure that spanning tree is connected, i.e., really a tree
        outputGateConnected = self.connectedStatic(self.n, var_span_edge)
        outputGate.appendToInput(outputGateConnected)

        # gates which give edges after removing the spanning tree
        remainingEdges = [[AndGate(self.id(), [self.var_edge(u, v), -var_span_edge(u, v)]) if u < v else None for v in self.V] for u in self.V]

        def var_remainingEdge(u, v):
            return remainingEdges[min(u, v)][max(u, v)]

        # gates which indicate whether it has degree 1 in the remaining graph, otherwise degree 0 or 2
        hasDegree1 = []
        for u in self.V:
            uDegree1 = OrGate(self.id(), [])
            for v in self.V:
                if v == u:
                    continue
                vIsOnlyNeighbor = AndGate(self.id(), [var_remainingEdge(u, v)] + [-var_remainingEdge(u, w) for w in self.V if w not in [u, v]])
                uDegree1.appendToInput(vIsOnlyNeighbor)
            hasDegree1.append(uDegree1)

        # ensure that only degree 1 vertices are adjacent to other degree 1 vertices
        for u, v in permutations(self.V, 2):
            # if u has degree 1 and there is a non-spanning tree edge to v, then also v has degree 1
            outputGate.appendToInput(OrGate(self.id(), [-hasDegree1[u], -var_remainingEdge(u, v), hasDegree1[v]]))
        return outputGate

    def treewidthVertexCritical(self, t, version=1):
        """Ensure that deleting any vertex results in treewidth t - 1"""
        outputGate = AndGate(self.id(), [])

        outputGate.appendToInput(self.treewidthAtMost(t, version))
        FalseGate = OrGate(self.id(), [])
        for v in self.V:

            def var_edge_adapted(i, j):
                if v in [i, j]:
                    return FalseGate
                else:
                    return self.var_edge(i, j)

            outputGate.appendToInput(self.treewidthAtMost(t - 1, version, var_edge=var_edge_adapted))
        return outputGate

    def treewidthCritical(self, t, version=1):
        """Ensure that every minor has treewidth t - 1"""
        outputGate = AndGate(self.id(), [])
        outputGate.appendToInput(self.treewidthAtMost(t, version))
        FalseGate = OrGate(self.id(), [])

        # ensure minimum degree 1 to avoid isolated vertices
        for u in self.V:
            outputGate.appendToInput(OrGate(self.id(), [self.var_edge(u, v) for v in self.V if v != u]))

        E = list(combinations(self.V, 2))  # [(0,i) for i in self.V if i > 0] #
        # deleting edge
        for u, v in E:
            # only applicable if edge is present, so if not present just make empty graph
            edgeTable = [[None for _ in self.V] for _ in self.V]

            for i, j in combinations(self.V, 2):
                if (i, j) == (u, v):
                    edgeTable[j][i] = edgeTable[i][j] = FalseGate
                else:
                    edgeTable[j][i] = edgeTable[i][j] = AndGate(
                        self.id(), [self.var_edge(i, j), self.var_edge(u, v)]
                    )  # empty graph if the edge is not present

            def var_edge_adapted(i, j):
                return edgeTable[i][j]

            outputGate.appendToInput(self.treewidthAtMost(t - 1, version, var_edge=var_edge_adapted))

        # merging edge
        for u, v in E:
            # add neighbor from v to u, and delete all neighbors from v

            edgeTable = [[None for _ in self.V] for _ in self.V]

            for i, j in combinations(self.V, 2):
                if v in [i, j]:
                    edgeTable[j][i] = edgeTable[i][j] = FalseGate
                elif u == i:
                    adjacent = OrGate(self.id(), [self.var_edge(i, j), self.var_edge(v, j)])
                    edgeTable[j][i] = edgeTable[i][j] = AndGate(self.id(), [adjacent, self.var_edge(u, v)])
                elif u == j:
                    adjacent = OrGate(self.id(), [self.var_edge(i, j), self.var_edge(i, v)])
                    edgeTable[j][i] = edgeTable[i][j] = AndGate(self.id(), [adjacent, self.var_edge(u, v)])
                else:
                    edgeTable[j][i] = edgeTable[i][j] = AndGate(
                        self.id(), [self.var_edge(i, j), self.var_edge(u, v)]
                    )  # empty graph if the edge is not present

            def var_edge_adapted(i, j):
                return edgeTable[i][j]

            outputGate.appendToInput(self.treewidthAtMost(t - 1, version, var_edge=var_edge_adapted))
        return outputGate

    def treewidthAtMost(self, t, version=2, var_edge=None):
        """
        Add an encoding which ensures that the treewidth is at most t from the graph given by var_edge.
        Returns a Gate capturing the constraint and list of the variables used.
        """
        outputGate = AndGate(self.id(), [])

        if var_edge == None:
            var_edge = self.var_edge
        ordering = [[self.id() if u < v else None for v in self.V] for u in self.V]  # if ture  then u < v else v < u, so automatically

        def ord(u, v):
            if u < v:
                return ordering[u][v]
            else:
                return -ordering[v][u]  # swap sign

        # transitivity
        for u, v, w in permutations(self.V, 3):
            outputGate.appendToInput(OrGate(self.id(), [-ord(u, v), -ord(v, w), ord(u, w)]))

        # TODO eventually symmetry breaking on twins, i.e., must be sorted by vertex labeling and can also be deleted right after each other

        if version == 1:  # more comapact but additional variables
            # arcs are the edges in the elemination orderings
            arcs = [[self.id() if u < v else None for v in self.V] for u in self.V]

            def var_arc(u, v):
                return arcs[min(u, v)][max(u, v)]

            for u, v in combinations(self.V, 2):
                # if edge then also arc
                outputGate.appendToInput(OrGate(self.id(), [-var_edge(u, v), var_arc(u, v)]))

                # closed
                for w in self.V:
                    if w in [u, v]:
                        continue
                    # if w adjacent to u and w and comes before both in the ordering then u,v must have an arc
                    outputGate.appendToInput(OrGate(self.id(), [-var_arc(w, u), -var_arc(w, v), -ord(w, u), -ord(w, v), var_arc(u, v)]))

            # for each vertex check degree to larger ones
            for u in self.V:
                postNeighbors = []
                for v in self.V:
                    if v == u:
                        continue
                    adjacentAndLater = AndGate(self.id(), [var_arc(u, v), ord(u, v)])
                    postNeighbors.append(adjacentAndLater)

                seqCounter(postNeighbors, countUpto=t, atMost=t, outputGate=outputGate, vPool=self)

        if version == 2:
            # second version: arc represent the smallest supergraph resulting from the elemination ordering.
            maxRounds = self.n - (t - 1)
            arcs_lvl = [
                [[OrGate(self.id(), []) if u < v else None for v in self.V] for u in self.V] for _ in range(maxRounds)
            ]  # graph at specific level
            arcs_lvl[0] = [[var_edge(u, v) if u < v else None for v in self.V] for u in self.V]  # first level are just the edges
            for r in range(1, maxRounds):
                for u, v in combinations(self.V, 2):
                    arcs_lvl[r][u][v].appendToInput(arcs_lvl[r - 1][u][v])  # edges remain edges
                    for w in self.V:
                        if w in [u, v]:
                            continue
                        # check whether w is a common neighbor in the previous graph
                        predecessorAndCommonNeighbor = AndGate(
                            self.id(), [ord(w, u), ord(w, v), arcs_lvl[r - 1][min(w, u)][max(w, u)], arcs_lvl[r - 1][min(w, v)][max(w, v)]]
                        )
                        arcs_lvl[r][u][v].appendToInput(predecessorAndCommonNeighbor)

            for u in self.V:
                postNeighbors = []
                for v in self.V:
                    if v == u:
                        continue
                    adjacentAndLater = AndGate(self.id(), [arcs_lvl[-1][min(u, v)][max(u, v)], ord(u, v)])
                    postNeighbors.append(adjacentAndLater)
                seqCounter(postNeighbors, countUpto=t, atMost=t, outputGate=outputGate, vPool=self)

        return outputGate

    def cardinialityConstraint(self, variables, countUpTo, counterType="sequential", atMost=None, atLeast=None):
        """Adds encoding to self (without using additional variables, but only gates) and returns the gates representing the counts"""
        if counterType == "sequential":
            seqCounter(variables, countUpTo, self, self.outputGate, atMost=atMost, atLeast=atLeast)
        else:
            print("Error: only support sequential counter")
            exit(1)

    def negateFormula(self):
        """
        Negate the gate assuming it being an AndGate
        """
        # self.outputGate = AndGate(self.id(), [NegatedGate(self.outputGate)])
        assert isinstance(self.outputGate, AndGate)
        self.outputGate = OrGate(self.id(), [NegatedGate(g) for g in self.outputGate.getInput()])

    def add_constraints_by_arguments(self, args):
        """Add constraints based on args given by default parser"""
        g = self
        # if g.DEBUG:
        #     print("Arguments:", args)
        #     stdout.flush()

        if args.chi_low:
            g = NegatedGate(self.colorable(args.chi_low - 1))
            self.addUniversalGate(g)

        if args.chi_upp:
            g = self.colorable(args.chi_upp)
            self.addExistentialGate(g)

        if args.Delta_upp:
            g = self.maxDegree(args.Delta_upp)
            self.addExistentialGate(g)
        if args.delta_low:
            g = self.minDegree(args.delta_low)
            self.addExistentialGate(g)

        if args.num_edges_upp:
            self.numEdgesUpp(args.num_edges_upp, args.counter)
        if args.num_edges_low:
            self.numEdgesLow(args.num_edges_low, args.counter)

        if args.Ck_free:
            g = self.ckFree(args.Ck_free)
            self.addExistentialGate(g)

        if args.tree_width_upp_version1:
            g = self.treewidthAtMost(args.tree_width_upp_version1, version=1)
            self.addExistentialGate(g)

        if args.tree_width_upp_version2:
            g = self.treewidthAtMost(args.tree_width_upp_version2, version=2)
            self.addExistentialGate(g)

        if args.tree_width_low_version1:
            g = self.treewidthAtMost(args.tree_width_low_version1, -1, version=1)
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.tree_width_low_version2:
            g = self.treewidthAtMost(args.tree_width_low_version2 - 1, version=2)
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.tree_width_vertex_critical:
            g = self.treewidthVertexCritical(args.tree_width_vertex_critical)
            self.addExistentialGate(g)

        if args.tree_width_critical:
            g = self.treewidthCritical(args.tree_width_critical)
            self.addExistentialGate(g)

        if args.non_010:
            g = self.colorable010()
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.min_domination_number:
            g = self.dominatingSet(args.min_domination_number - 1)
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.kochen_specker:
            """Ensure that graph is not 010 colorable"""
            # universal part
            triangleVars = {
                (u, v, w): AndGate(self.id(), [self.var_edge(u, v), self.var_edge(v, w), self.var_edge(u, w)]) for u, v, w in combinations(self.V, 3)
            }

            g = self.colorable010(triangleVars)
            g = NegatedGate(g)
            self.addUniversalGate(g)

            # existential part
            result = AndGate(self.id(), inputs=[])

            # each vertex in a triangle
            for u in self.V:
                g = OrGate(self.id(), [triangleVars[tuple(sorted((u, v, w)))] for v, w in combinations(self.V, 2) if u not in [v, w]])
                result.appendToInput(g)

            # max chromatic number 4
            g = self.colorable(4)
            result.appendToInput(g)

            # min degree 3
            g = self.minDegree(3)
            result.appendToInput(g)

            # no 4 cycle
            g = self.ckFree(4)
            result.appendToInput(g)

            self.addExistentialGate(result)

        if args.no_subsuming_neighborhoods:
            g = self.noSubsumingNeighboorhoods()
            self.addExistentialGate(g)

        if args.mtf:
            g = self.mtf()
            self.addExistentialGate(g)

        if args.non_3_edge_colorable:
            g = self.threeEdgeColorable()
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.cubic:
            g = self.cubic()
            self.addExistentialGate(g)

        if args.two_connected:
            g = self.twoConnected()
            self.addExistentialGate(g)

        if args.folkmann:
            g = self.triangleEdgeColoring()
            g = NegatedGate(g)
            self.addUniversalGate(g)

            g = self.maxCliqueSize(args.folkmann - 1)
            self.addExistentialGate(g)

            g = self.vertexInTwoTriangles()
            self.addExistentialGate(g)

            if args.folkmann == 4:
                g = self.maximal4cliqueFree()
                self.addExistentialGate(g)

        if args.maximal_4_clique_free:
            g = self.maximal4cliqueFree()
            self.addExistentialGate(g)

        if args.max_independent_set:
            g = self.maxIndependentSet(args.max_independent_set)
            self.addExistentialGate(g)

        if args.three_decomposition_conjecture:
            g = self.specialSpanningTree()
            g = NegatedGate(g)
            self.addUniversalGate(g)

            g = self.cubic()
            self.addExistentialGate(g)

        if args.bipartite:
            g = self.bipartite()
            self.addExistentialGate(g)

        if args.min_girth_compact:
            self.minGirthCompact(args.min_girth_compact)

        if args.domination_conjecture:
            import math

            minDominationNumber = math.ceil(self.n / 3) + 1
            g = self.dominatingSet(minDominationNumber - 1)
            g = NegatedGate(g)
            self.addUniversalGate(g)

            g = self.cubic()
            self.addExistentialGate(g)

        if args.connected_static:
            g = self.connectedStatic(self.n)
            self.addExistentialGate(g)

        if args.canonical_qbf:
            g = self.lexSmallerMatrix(useColexOrdering=False)
            g = NegatedGate(g)
            self.addUniversalGate(g)

        if args.canonical_qbf_colex:
            g = self.lexSmallerMatrix(useColexOrdering=True)
            g = NegatedGate(g)
            self.addUniversalGate(g)

    # ------------degree encodings--------------

    def minDegree(self, delta):
        """Minimum degree at least delta"""
        outputGate = AndGate(self.id(), [])
        g = self
        for u in g.V:
            seqCounter([g.var_edge(u, v) for v in g.V if v != u], delta, self, outputGate, atLeast=delta)
        return outputGate

    def maxDegree(self, delta):
        """Maximum degree at most delta"""
        outputGate = AndGate(self.id(), [])
        g = self
        for u in g.V:
            seqCounter([g.var_edge(u, v) for v in g.V if v != u], delta, self, outputGate, atMost=delta)
        return outputGate

    # -------------------number of edges ------------------

    def numEdgesUpp(self, m):
        """Upperbound on edges"""
        g = self.outputGate
        seqCounter([g.var_edge(u, v) for u, v in combinations(self.V, 2)], m, self, g, atMost=m)

    def numEdgesLow(self, m):
        """Lowerbound on edges"""
        g = self.outputGate
        seqCounter([g.var_edge(u, v) for u, v in combinations(self.V, 2)], m, self, g, atLeast=m)

    def ckFree(self, k):
        """Forbid k-cycles (C_k) as subgraphs"""
        outputGate = AndGate(self.id(), [])
        for cycle in permutations(self.V, k):
            if cycle[0] != min(cycle):
                continue
            if cycle[1] > cycle[-1]:
                continue
            outputGate.appendToInput(
                OrGate(self.id(), [-self.var_edge(cycle[i], cycle[(i + 1) % k]) for i in range(k)])
            )  # at least one edge absent from potential cycle
        return outputGate

    def noSubsumingNeighboorhoods(self):
        outputGate = AndGate(self.id(), [])
        for u, v in permutations(self.V, 2):
            different = OrGate(self.id(), [self.var_edge(u, v)])
            for w in self.V:
                if w in [u, v]:
                    continue
                different.appendToInput(AndGate(self.id(), [+self.var_edge(u, w), -self.var_edge(v, w)]))
            outputGate.appendToInput(different)
        return outputGate

    def mtf(self):
        """Ensure that maximal triangle free i.e., no 3 cycle and inserting any additional edge results in a triangle"""
        outputGate = AndGate(self.id(), [])
        for u, v, w in combinations(self.V, 3):
            outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v), -self.var_edge(v, w), -self.var_edge(u, w)]))
        for u, v in combinations(self.V, 2):
            # edge not present implies that the have a common neighbor
            outputGate.appendToInput(
                OrGate(
                    self.id(),
                    [self.var_edge(u, v)] + [AndGate(self.id(), [self.var_edge(u, w), self.var_edge(v, w)]) for w in self.V if w not in [u, v]],
                )
            )
        return outputGate

    def greedyColoring(self):
        """Ensure that it can be colored by greedy coloring"""

        # Create variables indicating ordering of vertices
        ordering = [
            [self.id() if u < v else None for v in self.V] for u in self.V
        ]  # if true  then u < v else v < u, so automatically totally ordered
        outputGate = AndGate(self.id(), [])

        def isSmaller(u, v):
            if u < v:
                return ordering[u][v]
            else:
                return -ordering[v][u]

        # compute the color for each vertex based on the ordering (Difficult to represent in this format)
        # More precisely color can not be a gate, because dependent from each other based on ordering

        exit(1)  # not implemented yet
        return outputGate

    def threeEdgeColorable(self):
        """Ensure that edges can be colored by three colors such that incident edges don't have the same color"""
        C = range(3)
        coloring = [[[self.id() if u < v else None for c in C] for v in self.V] for u in self.V]
        outputGate = AndGate(self.id(), [])
        # each edge has at least one color
        for u, v in combinations(self.V, 2):
            outputGate.appendToInput(OrGate(self.id(), coloring[u][v]))

        # no two incident edges have the same color
        for u1, v1 in combinations(self.V, 2):
            for u2, v2 in combinations(self.V, 2):
                if u1 > u2 or (u1 == u2 and v1 >= v2):
                    continue  # only in one direction
                if u1 not in [u2, v2] and v1 not in [u2, v2]:
                    continue  # must be incident

                # print(u1, v1, u2, v2)

                # if both edges present than not both have the same color
                for c in C:
                    outputGate.appendToInput(
                        OrGate(self.id(), [-self.var_edge(u1, v1), -self.var_edge(u2, v2), -coloring[u1][v1][c], -coloring[u2][v2][c]])
                    )
        return outputGate

    def cubic(self):
        """Ensure that graph is cubic"""
        outputGate = AndGate(self.id(), [])
        degree = 3
        for u in self.V:
            seqCounter([self.var_edge(u, v) for v in self.V if v != u], degree, self, outputGate, atMost=degree, atLeast=degree)
        return outputGate

    def twoConnected(self):
        """Deleting any vertex doesn't make the graph disconnected"""
        outputGate = AndGate(self.id(), [])
        for u in self.V:

            def m(v):  # remapping of the remaining vertices
                return v if v < u else v + 1

            def var_edge_remaining(v1, v2):
                return self.var_edge(m(v1), m(v2))

            outputGateU = self.connectedStatic(self.n - 1, var_edge_remaining)
            outputGate.appendToInput(outputGateU)

        return outputGate

    def triangleEdgeColoring(self):
        """Edges can be colored by two colors such that no monochromatic triangle"""
        outputGate = AndGate(self.id(), [])
        # create ids for each edge representing the color
        coloring = [[self.id() if u < v else None for v in self.V] for u in self.V]
        # no monochromatic triangle
        for u, v, w in combinations(self.V, 3):
            outputGate.appendToInput(
                OrGate(
                    self.id(), [-self.var_edge(u, v), -self.var_edge(u, w), -self.var_edge(v, w), -coloring[u][v], -coloring[v][w], -coloring[u][w]]
                )
            )
            outputGate.appendToInput(
                OrGate(self.id(), [-self.var_edge(u, v), -self.var_edge(u, w), -self.var_edge(v, w), coloring[u][v], coloring[v][w], coloring[u][w]])
            )
        return outputGate

    def maxCliqueSize(self, k):
        """Maximum clique size is k"""
        outputGate = AndGate(self.id(), [])
        for clique in combinations(self.V, k + 1):  # forbid k + 1 clique
            outputGate.appendToInput(OrGate(self.id(), [-self.var_edge(u, v) for u, v in combinations(clique, 2)]))
        return outputGate

    def maxIndependentSet(self, k):
        """Maximum size of an independent set is k"""
        outputGate = AndGate(self.id(), [])
        for clique in combinations(self.V, k + 1):  # forbid k + 1 independent set
            outputGate.appendToInput(OrGate(self.id(), [self.var_edge(u, v) for u, v in combinations(clique, 2)]))
        return outputGate

    def vertexInTwoTriangles(self):
        """Each vertex is in two triangles sharing an edge"""
        outputGate = AndGate(self.id(), [])

        for u in self.V:
            structureIsPresent = OrGate(self.id(), [])
            for v in self.V:  # the second vertex in the shared edge
                if v == u:
                    continue
                for w1, w2 in combinations(set(self.V) - {v, u}, 2):
                    structureIsPresent.appendToInput(
                        AndGate(
                            self.id(), [self.var_edge(u, v), self.var_edge(u, w1), self.var_edge(u, w2), self.var_edge(v, w1), self.var_edge(v, w2)]
                        )
                    )

            outputGate.appendToInput(structureIsPresent)

        return outputGate

    def maximal4cliqueFree(self):
        """Inserting any additional edge results in a 4-clique"""
        outputGate = self.maxCliqueSize(3)

        for u, v in combinations(self.V, 2):
            checkGate = OrGate(self.id(), [self.var_edge(u, v)])
            for w1, w2 in combinations(set(self.V) - {u, v}, 2):
                checkGate.appendToInput(
                    AndGate(self.id(), [self.var_edge(x1, x2) for x1, x2 in combinations([u, v, w1, w2], 2) if (x1, x2) != (u, v)])
                )  # almost clique, so reason to not add edge
            outputGate.appendToInput(checkGate)
        return outputGate

    def solve(self, allGraphs=False, hideGraphs=False, qcirFile=None, args_SMS="", forwarding_args=[], graph6_format=False) -> None:
        """Solve the formula, given the encoding, using SMS.

        :param allGraphs: Enumerate all satisfying graphs. Default value = False
        :param hideGraphs: Count all satisfying graphs. Default value = False
        :param qcirFile: Write constraints here, use a temporary if None. Default value = None
        :param forwarding_args: Forward these arguments to SMS. Default value = []
        :param graph_format: Toggle output format. Choices are [graph6, edge_list]. Default value = edge-list
        :param args_SMS: Forward these arguments to SMS. Deprecated, use `forwarding_args` instead. Default value = ""

        """
        if qcirFile == None:
            import time

            qcirFile = f"./temp{os.getpid()}_t{time.time()}.enc"  # TODO use tempfile module
        with open(qcirFile, "w") as f:
            self.print_qcir(f)  # write script to temporary file

        program = "smsd" if self.directed else "smsg"  # we expect these binaries to be on PATH

        # add arguments

        if allGraphs:
            self.paramsSMS["all-graphs"] = ""
        if hideGraphs:
            self.paramsSMS["hide-graphs"] = ""

        python_args_SMS = " ".join(f"--{param} {value}" for param, value in self.paramsSMS.items())

        sms_command = "time " if self.DEBUG else ""
        sms_command += (
            f"{program} {python_args_SMS} {args_SMS} --qcir-file {qcirFile}"  # TODO eventually parse args_SMS to allow to override deault arguments
        )
        for arg in forwarding_args:
            sms_command += f" '{arg}'"

        if self.DEBUG:
            print("running the command: ", sms_command)

        if graph6_format:
            import networkx as nx
            from ast import literal_eval

            for line in os.popen(sms_command).read().split("\n"):
                if line and line[0] == "[":
                    edges = literal_eval(line)
                    print(nx.to_graph6_bytes(nx.Graph(edges), header=False).decode(), end="")
                elif self.DEBUG:
                    print(line, end="\n")
        else:
            os.system(sms_command)

        os.system(f"rm {qcirFile}")  # cleanup

    def solveArgs(self, args, forwarding_args) -> None:
        """Wrapper for solving using arguments provided by argsParser

        :param args: Arguments for PySMS
        :param forwarding_args: Arguments to be forwarded to SMS

        """
        self.solve(
            allGraphs=args.all_graphs,
            hideGraphs=args.hide_graphs,
            qcirFile=args.qcir_file,
            args_SMS=args.args_SMS,
            forwarding_args=forwarding_args,
            graph6_format=args.graph6_format,
        )


if __name__ == "__main__":
    args, forwarding_args = getDefaultParser().parse_known_args()
    if forwarding_args:
        print("WARNING: Unknown arguments for python script which are forwarded to SMS:", forwarding_args, file=stderr)
    b = GraphEncodingBuilder(args.vertices, directed=args.directed)
    b.add_constraints_by_arguments(args)
    if args.print_cnf:
        with open(args.print_cnf, "w") as cnf_fh:
            b.print_dimacs(cnf_fh)
    if args.print_qcir:
        with open(args.print_qcir, "w") as qcir_fh:
            b.print_qcir(qcir_fh)
    else:
        b.solveArgs(args, forwarding_args)

    # -------------------test cases-------------------
    """
    b = GraphEncodingBuilder(10)
    g1 = AndGate(200, [1, -2])
    g2 = OrGate(201, [1, -2, b.var_edge(1, 4)])
    b.outputGate.inputs = [g1, g2]
    b.print_qcir()

    print("----------------")
    b = GraphEncodingBuilder(5)
    b.colorable(3)
    b.negateFormula()
    b.print_qcir()

    print("----------------")
    b = GraphEncodingBuilder(5)
    b.treewidthAtMost(3)
    b.print_qcir()

    # -------------simple cardinality constraints sanity checks
    testCard = 0
    if testCard == 1:
        b = GraphEncodingBuilder(5)
        b.cardinialityConstraint([b.var_edge(0, i) for i in [1, 2, 3, 4]], 2, atMost=2)
        # b.outputGate.appendToInput(AndGate(b.id(), [b.var_edge(0, 1), b.var_edge(0, 2), b.var_edge(0, 3)]))
        b.print_qcir()
    if testCard == 2:
        b = GraphEncodingBuilder(5)
        b.cardinialityConstraint([b.var_edge(0, i) for i in [1, 2, 3, 4]], 2, atMost=2)
        b.outputGate.appendToInput(AndGate(b.id(), [b.var_edge(0, 2), b.var_edge(0, 3)]))
        b.print_qcir()
    if testCard == 3:
        b = GraphEncodingBuilder(5)
        b.cardinialityConstraint([b.var_edge(0, i) for i in [1, 2, 3, 4]], 2, atMost=2)
        b.outputGate.appendToInput(AndGate(b.id(), [b.var_edge(0, 1), b.var_edge(0, 2), b.var_edge(0, 3)]))
        b.print_qcir()
    if testCard == 4:
        b = GraphEncodingBuilder(5)
        b.cardinialityConstraint([b.var_edge(0, i) for i in [1, 2, 3, 4]], 2, atMost=2)
        b.outputGate.appendToInput(AndGate(b.id(), [b.var_edge(0, 2), b.var_edge(0, 3), b.var_edge(0, 4)]))
        b.print_qcir()

    # --------------------test regarding tree width------------------------
    testTreeWidth = 0

    # K_4 has treewidth at most 3
    if testTreeWidth == 1:
        b = GraphEncodingBuilder(4)
        # make clique
        for u, v in combinations(b.V, 2):
            b.outputGate.appendToInput(b.var_edge(u, v))
        b.treewidthAtMost(3)
        b.print_qcir()

    # K_4 doesn't have treewidth at most 2
    if testTreeWidth == 2:
        b = GraphEncodingBuilder(4)
        # make clique
        for u, v in combinations(b.V, 2):
            b.outputGate.appendToInput(b.var_edge(u, v))
        b.treewidthAtMost(2)
        b.print_qcir()

    # K_4 minus one edge has treewidth at most 2
    if testTreeWidth == 3:
        b = GraphEncodingBuilder(4)
        # make clique
        for u, v in combinations(b.V, 2):
            if (u, v) == (0, 1):
                continue
            b.outputGate.appendToInput(b.var_edge(u, v))
        b.treewidthAtMost(2)
        b.print_qcir()

    if testTreeWidth == 4:
        b = GraphEncodingBuilder(10)
        # make clique
        for u, v in combinations(b.V, 2):
            b.outputGate.appendToInput(b.var_edge(u, v))
        b.treewidthAtMost(9)
        b.print_qcir()

    if False:
        b = GraphEncodingBuilder(10)
        # make clique
        for u, v in combinations(b.V, 2):
            b.outputGate.appendToInput(b.var_edge(u, v))
        b.treewidthAtMost(9)
        b.print_dimacs()

    args = getDefaultParser().parse_args()
    b = GraphEncodingBuilder(args.vertices)
    b.add_constraints_by_arguments(args)

    if args.negate_output_gate:
        b.negateFormula()
    if args.print_dimacs:
        b.print_dimacs()
    else:
        b.print_qcir()
    """
