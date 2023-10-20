""" 
This file can either directly be used for generating all graphs described by some predefined properties which can be selected with flags,
or it can be use as starting point for constructing an encoding for graphs on your own.
"""

from itertools import *
from pysms.counters import *
import argparse
import os
from sys import *
from ast import literal_eval


def getDefaultParser():
    """ """
    parser = argparse.ArgumentParser()  # default parser for some properties

    main_args = parser.add_argument_group(title="Main arguments", description="The number of vertices is mandatory, everything else is optional")
    solve_args = parser.add_argument_group(title="Solver options")
    constraint_args = parser.add_argument_group(title="Graph constraints", description="A set of pre-defined constraints for common applications, including applications from SMS papers")

    main_args.add_argument("--vertices", "-v", type=int, required=True, help="number of vertices")
    main_args.add_argument("--cnf-file", type=str, help="store the generated encoding here")
    main_args.add_argument("--directed", "->", action="store_true", help="search for directed graphs")
    main_args.add_argument("--multigraph", "-m", type=int, help="search for a multigraph")
    main_args.add_argument("--underlying-graph", action="store_true", help="consider the underlying undirected graph for directed graphs")
    main_args.add_argument("--static-partition", action="store_true", help="specify a statically enforced partial vertex ordering (respected by SMS)")
    main_args.add_argument("--counter", choices=["sequential", "totalizer"], default="sequential", help="the CNF encoding for cardinality constraints")
    main_args.add_argument("--DEBUG", type=int, default=1, help="debug level")

    # solve options
    solve_args.add_argument("--no-solve", action="store_true", help="don't run SMS, only generate the constraints and output either to stdout, or to a file specified by --cnf-file")
    solve_args.add_argument("--all-graphs", "-a", action="store_true", help="generate all graphs (without this, the solver will exit after the first solution)")
    solve_args.add_argument("--hide-graphs", "-hg", action="store_true", help="do not display graphs (meant as a counting functionality, though the graphs still need to be enumerated)")
    solve_args.add_argument("--args-SMS", type=str, default="", help="command line to be appended to the call to smsg/smsd (see src/main.cpp or README.md)")
    solve_args.add_argument("--graph-format", choices=['edge-list','graph6'], default=None, help="output format of graphs")

    # number of edges
    constraint_args.add_argument("--num-edges-upp", "-E", type=int, help="upper bound on the maximum number of edges")
    constraint_args.add_argument("--num-edges-low", "-e", type=int, help="lower bound on the minimum number of edges")

    # degree restrictions
    constraint_args.add_argument("--Delta-upp", "-D", type=int, help="upper bound on the maximum degree")
    constraint_args.add_argument("--delta-low", "-d", type=int, help="lower bound on the minimum degree")
    constraint_args.add_argument("--even-degrees", action="store_true", help="all degrees should be even")
    constraint_args.add_argument("--no-subsuming-neighborhoods", action="store_true", help="ensure that N(v) ⊈ N(u) for any vertex pair u,v")
    constraint_args.add_argument("--degree-partition", action="store_true", help="sort vertices by degree and only apply SMS on vertices with same degree")

    # chromatic number
    constraint_args.add_argument("--chi-upp", type=int, help="upper bound on the chromatic number")
    constraint_args.add_argument("--chi-low", type=int, help="lower bound on the chromatic number (not encoded to CNF, needs SMS)")

    # cycles
    constraint_args.add_argument("--Ck-free", type=int, help="forbid the k-cycle C_k as (non-induced) subgraph")
    constraint_args.add_argument("--mtf", action="store_true", help="search for maximal triangle-free graphs (where adding any edge creates a triangle)")
    constraint_args.add_argument("--girth", type=int, help="lower bound on girth")
    constraint_args.add_argument("--girth-compact", type=int, help="lower bound on girth, more compact encoding")

    # ramsey theory
    constraint_args.add_argument("--alpha-upp", type=int, help="maximum size of an independent set")
    constraint_args.add_argument("--omega-upp", type=int, help="maximum size of a clique")
    constraint_args.add_argument("--ramsey", nargs=2, type=int, help="(a, w) means no independent set of size a and no clique of size w")

    # other
    constraint_args.add_argument("--planar-kuratowski", "--planar", "-p", action="store_true", help="generate only planar graphs (not encoded to CNF, needs SMS)")
    constraint_args.add_argument("--connectivity-low", "-c", "--kappa-low", default=0, type=int, help="lower bound on vertex connectivity")  # TODO handle in SMS
    constraint_args.add_argument("--diam2-critical", action="store_true", help="assert a diameter-2-critical graph")

    return parser


class IDPool:
    """A class for returning the next free id starting with 1"""

    def __init__(self, start_from=1) -> None:
        self.nextId = start_from

    def id(self):
        """Returns the next free variable"""
        x = self.nextId
        self.nextId += 1
        return x


def CNF_OR(ins, out) -> list[int]:
    """Returns a list of clauses ensuring that `out` is equivalent to the disjunction of the literals in `ins`

    :param ins: input variable ids
    :param out: a (fresh) variable id

    """
    return [[-out] + ins] + [[out, -x] for x in ins]


def CNF_AND(ins, out) -> list[int]:
    """Returns a list of clauses ensuring that `out` is equivalent to the conjunction of the literals in `ins`

    :param ins: input variable ids
    :param out: a (fresh) variable id

    """
    return [[out] + [-x for x in ins]] + [[-out, x] for x in ins]


DEFAULT_COUNTER = "sequential"


class GraphEncodingBuilder(IDPool, list):
    """A class for building an encoding in the context of SMS.
    The IDPool gives the next free variable whilst the list contains the clauses


    """

    def __init__(self, n, directed=False, multiGraph=None, staticInitialPartition=False, underlyingGraph=False, DEBUG=0):
        super().__init__()
        self.directed = directed
        self.V = list(range(n))
        self.n = n
        self.DEBUG = DEBUG
        self.varStaticInitialPartition = None

        self.paramsSMS = {"vertices": self.n, "print-stats": True, "frequency": 30}  # default params

        # order in which variables are assigned is import !!!
        if directed:
            self.varEdgeDirectedTable = [[None for _ in self.V] for _ in self.V]
            for v, u in permutations(self.V, 2):
                self.varEdgeDirectedTable[v][u] = self.id()
        elif multiGraph:
            self.paramsSMS["multi-graph"] = multiGraph
            self.varEdgeMultiTable = [[[None for _ in self.V] for _ in self.V] for _ in range(multiGraph)]
            for i in range(multiGraph):
                for v, u in combinations(self.V, 2):
                    self.varEdgeMultiTable[i][v][u] = self.varEdgeMultiTable[i][u][v] = self.id()
            self.varEdgeTable = self.varEdgeMultiTable[0]  # allows arguing over the first graph
        else:
            self.varEdgeTable = [[None for _ in self.V] for _ in self.V]
            for v, u in combinations(self.V, 2):
                self.varEdgeTable[v][u] = self.varEdgeTable[u][v] = self.id()

        if staticInitialPartition:
            self.varStaticInitialPartition = [[None for _ in self.V] for _ in self.V]
            for v, u in combinations(self.V, 2):
                self.varStaticInitialPartition[v][u] = self.varStaticInitialPartition[u][v] = self.id()

            for u, v, w in combinations(self.V, 3):
                # ensure that partition is well defined, i.e., elements in the middle are equal
                self.append([-self.var_partition(u, w), self.var_partition(u, v)])
                self.append([-self.var_partition(u, w), self.var_partition(v, w)])
                # transitive
                self.append([-self.var_partition(u, v), -self.var_partition(v, w), self.var_partition(u, w)])

        if underlyingGraph:
            self.varEdgeTable = [[None for _ in self.V] for _ in self.V]
            for v, u in combinations(self.V, 2):
                self.varEdgeTable[v][u] = self.varEdgeTable[u][v] = self.id()
                self.CNF_OR_APPEND(
                    [self.varEdgeDirectedTable[v][u], self.varEdgeDirectedTable[u][v]], self.varEdgeTable[v][u]
                )  # relelation between edge variables and directed edge variables, i.e., undirected is the underlying graph

    def var_edge(self, u, v) -> int:
        """Get the propositional variable associated with the undirected edge {u,v}

        :param u: a vertex
        :param v: a different vertex
        :returns: the undirected edge variable e_{u,v}

        """
        return self.varEdgeTable[u][v]

    def var_edge_multi(self, i, u, v) -> int:
        """Get the propositional variable associated with the undirected edge {u,v} of the i-th graph

        :param i: the graph index
        :param v: a vertex
        :param u: a different vertex
        :returns: the undirected multigraph edge variable e^i_{u,v}

        """
        return self.varEdgeMultiTable[i][u][v]

    def var_edge_dir(self, u, v) -> int:
        """Get the propositional variable associated with the directed edge (u,v)

        :param u: a vertex
        :param v: a different vertex
        :returns: the directed edge variable e_{u,v}

        """
        return self.varEdgeDirectedTable[u][v]

    def var_partition(self, u, v) -> int:
        """Get the propositional variable which holds whether u and v are in the same partition

        :param u: a vertex
        :param v: a different vertex
        :returns: the partition description variable

        """
        return self.varStaticInitialPartition[u][v]

    def solve(self, allGraphs=False, hideGraphs=False, cnfFile=None, args_SMS="", forwarding_args=[], graph_format=None) -> None:
        """Solve the formula, given the encoding, using SMS.

        :param allGraphs: Enumerate all satisfying graphs. Default value = False
        :param hideGraphs: Count all satisfying graphs. Default value = False
        :param cnfFile: Write constraints here, use a temporary if None. Default value = None
        :param forwarding_args: Forward these arguments to SMS. Default value = []
        :param graph_format: Toggle output format. Choices are [graph6, edge_list]. Default value = edge-list
        :param args_SMS: Forward these arguments to SMS. Deprecated, use `forwarding_args` instead. Default value = ""

        """
        if cnfFile == None:
            cnfFile = f"./temp{os.getpid()}.enc"  # TODO use tempfile module
        with open(cnfFile, "w") as cnf_fh:
            self.print_dimacs(cnf_fh)  # write script to temporary file

        program = "smsd" if self.directed else "smsg"  # we expect these binaries to be on PATH

        # add arguments

        if allGraphs:
            self.paramsSMS["all-graphs"] = ""
        if hideGraphs:
            self.paramsSMS["hide-graphs"] = ""
        if self.varStaticInitialPartition:
            self.paramsSMS["combine-static-dynamic"] = ""

        python_args_SMS = " ".join(f"--{param} {value}" for param, value in self.paramsSMS.items())

        sms_command = "time " if self.DEBUG else ""
        sms_command += f"{program} {python_args_SMS} {args_SMS} --dimacs {cnfFile}"  # TODO eventually parse args_SMS to allow to override
        for arg in forwarding_args:
            sms_command += f" '{arg}'"

        if self.DEBUG:
            print("running the command: ", sms_command)

        if graph_format:
            assert(graph_format in ['edge-list','graph6'])
            if graph_format == 'graph6':
                import networkx as nx
            for line in os.popen(sms_command).read().split("\n"):
                if line and line[0] == '[':
                    edges = literal_eval(line)
                    if graph_format == 'graph6':
                        print(nx.to_graph6_bytes(nx.Graph(edges),header=False).decode(),end="")
                    if graph_format == 'edge-list':
                        print(edges)
                elif self.DEBUG:
                    print(line,end="\n")
        else:
            os.system(sms_command)

        os.system(f"rm {cnfFile}") # cleanup


    def solveArgs(self, args, forwarding_args) -> None:
        """Wrapper for solving using arguments provided by argsParser

        :param args: Arguments for PySMS
        :param forwarding_args: Arguments to be forwarded to SMS

        """
        self.solve(allGraphs=args.all_graphs, hideGraphs=args.hide_graphs, cnfFile=args.cnf_file, args_SMS=args.args_SMS, forwarding_args=forwarding_args, graph_format=args.graph_format)

    # ------------------some utilies--------------------------

    def CNF_OR_APPEND(self, ins, out) -> None:
        """extends with an encoding of an or gate

        :param ins: input literal ids
        :param out: output gate variable id

        """
        self.extend(CNF_OR(ins, out))

    def CNF_AND_APPEND(self, ins, out) -> None:
        """extends with an encoding of an and gate

        :param ins: input literal ids
        :param out: output gate variable id

        """
        self.extend(CNF_AND(ins, out))

    def CNF_OR(self, ins) -> int:
        """returns a new variable which is true iff at least one of the literals in 'ins' is true

        :param ins: input literal ids
        :returns: a fresh Boolean variable equivalent to the disjunction of ins

        """
        out = self.id()
        self.extend(CNF_OR(ins, out))
        return out

    def CNF_AND(self, ins) -> int:
        """returns a new variable which is true iff all literals in 'ins' are true

        :param ins: input literal ids
        :returns: a fresh Boolean variable equivalent to the disjunction of ins

        """
        out = self.id()
        self.extend(CNF_AND(ins, out))
        return out

    def print_dimacs(self, outstream) -> None:
        """Print the current encoding to the given output stream

        :param outstream: print to here

        """
        print(f"p cnf {len(self)} {self.nextId-1}", file=outstream)
        for c in self:
            print(" ".join(str(x) for x in c), 0, file=outstream)

    def counterFunction(self, literals, countUpto, atMost=None, atLeast=None, counterType="sequential") -> list[int]:
        """A wrapper for `pysms.counters.counterFunction`: constraints are added to and the ids are provided by the builder object (self).

        :param literals: create a counter over these literals
        :param countUpto: extend the counter to this value
        :param atMost: enforce a cardinality upper bound. Default value = None
        :param atLeast: enforce a cardinality lower bound. Default value = None
        :param counterType: specify which counter encoding to use. Default value = "sequential"
        :returns: auxiliary counter variables

        """
        return counterFunction(literals, countUpto, self, self, atMost=atMost, atLeast=atLeast, type=counterType)

    # -------------------------ENCODINGS----------------------------------------

    def add_constraints_by_arguments(self, args) -> None:
        """Add constraints based on args given by default parser

        :param args: arguments parsed from the command line

        """
        g = self
        if g.DEBUG:
            print("Arguments:", args, file=stderr)
            stderr.flush()
        self.degreeBounds(self.V, args.delta_low, args.Delta_upp, encoding=args.counter)
        self.edgeBounds(self.V, args.num_edges_low, args.num_edges_upp, encoding=args.counter)
        if args.diam2_critical:
            self.diameter2critical()

        if args.Ck_free:
            self.ckFree(args.Ck_free)
        if args.girth:
            self.minGirth(args.girth)
        if args.girth_compact:
            self.minGirthCompact(args.girth_compact)

        if args.alpha_upp:
            self.maxIndependentSet(args.alpha_upp)
        if args.omega_upp:
            self.maxClique(args.omega_upp)
        if args.ramsey:
            self.maxIndependentSet(args.ramsey[0] - 1)
            self.maxClique(args.ramsey[1] - 1)

        if args.mtf:
            self.mtf()

        if args.no_subsuming_neighborhoods:
            self.noSubsumingNeighborhoods()
        if args.degree_partition:
            self.sort_vertices_by_degree()

        if args.chi_upp:
            self.maxChromaticNumber(args.chi_upp)

        if args.chi_low:
            self.paramsSMS["min-chromatic-number"] = args.chi_low

        if args.connectivity_low:
            self.minConnectivity(args.connectivity_low)

        if args.planar_kuratowski:
            self.paramsSMS["planar"] = 5  # DEFAULT planarity frequency

        if args.even_degrees:
            for u in self.V:
                shouldBe([+self.var_edge(u, v) for v in self.V if v != u], [i for i in self.V if i % 2 == 0], self, self, type=DEFAULT_COUNTER)

    # ------------degree encodings--------------

    def minDegree(self, delta, countertype=DEFAULT_COUNTER) -> list[list[int]]:
        """Minimum degree at least delta

        :param delta: degree lower bound
        :param countertype: specify a cardinality encoding. Default value = DEFAULT_COUNTER

        """
        return degreeBounds(self.V, delta, None, encoding=countertype)

    def maxDegree(self, delta, countertype=DEFAULT_COUNTER) -> list[list[int]]:
        """Maximum degree at most delta

        :param delta: degree upper bound
        :param countertype: specify a cardinality encoding. Default value = DEFAULT_COUNTER

        """
        return degreeBounds(self.V, None, delta, encoding=countertype)

    def degreeBounds(self, verts, lower, upper, within=False, encoding=DEFAULT_COUNTER) -> list[list[int]]:
        """Enforce that each of verts has degree between lower and upper (either of which can be None)

        :param verts: the set of vertices on which to operate
        :param lower: a lower bound for the minimum degree
        :param upper: an upper bound for the maximum degree
        :param within: whether to count the degree only within verts, or to the whole graph (Default value = False)
        :param encoding: the cardinality encoding to be used (Default value = DEFAULT_COUNTER)
        :returns: A list of lists of auxiliary variables, which say that the `i`-th vertex has degree at least `j`,
        or `None` if the bounds are meaningless.

        """
        if lower == None and upper == None:
            return None
        if lower and upper and lower > upper:
            self.append([])
            return None
        aux = []
        target_set = verts if within else self.V
        count_limit = upper if upper else lower
        for u in verts:
            aux.append(self.counterFunction([self.var_edge(u, v) for v in target_set if v != u], count_limit, upper, lower, counterType=encoding))
        return aux

    # -------------------number of edges ------------------

    def numEdgesUpp(self, m, countertype=DEFAULT_COUNTER) -> list[int]:
        """An upper bound on the number of edges

        :param m: an upper bound on the number of edges
        :param countertype: the cardinality encoding to be used (Default value = DEFAULT_COUNTER)

        """
        return self.edgeBounds(self.V, None, m, encoding=countertype)

    def numEdgesLow(self, m, countertype=DEFAULT_COUNTER) -> list[int]:
        """A lower bound on the number of edges

        :param m: a lower bound on the number of edges
        :param countertype: the cardinality encoding to be used (Default value = DEFAULT_COUNTER)

        """
        return self.edgeBounds(self.V, m, None, encoding=countertype)

    def edgeBounds(self, verts, lower, upper, encoding=DEFAULT_COUNTER) -> list[int]:
        """Enforce that each of verts has degree between lower and upper (either of which can be None)

        :param verts: the set of vertices between which to enforce edge bounds
        :param lower: lower bound on the number of edges
        :param upper: upper bound on the number of edges
        :param encoding: the cardinality encoding to be used (default DEFAULT_COUNTER)
        :returns: A list of auxiliary variables, which say that the number of edges is at least i,
        or None if the bounds are meaningless.

        """
        if lower == None and upper == None:
            return None
        if lower and upper and lower > upper and verts:
            self.append([])
            return None
        count_limit = upper if upper else lower
        return self.counterFunction([self.var_edge(u, v) for u, v in combinations(verts, 2)], count_limit, upper, lower, counterType=encoding)

    # ------------------------------------------------------

    def mtf(self) -> None:
        """Make the graph maximal triangle-free (adding any further edges will create triangles)"""
        self.ckFree(3)  # forbid triangles
        commonNeighbor = {(i, j, k): self.id() for i, j in combinations(self.V, 2) for k in self.V if k not in [i, j]}
        for i, j in combinations(self.V, 2):
            for k in self.V:
                if k in [i, j]:
                    continue
                self.CNF_AND_APPEND([self.var_edge(i, k), self.var_edge(j, k)], commonNeighbor[(i, j, k)])
            self.append([self.var_edge(i, j)] + [commonNeighbor[(i, j, k)] for k in self.V if k not in [i, j]])

    def noSubsumingNeighborhoods(self) -> None:
        """No neighborhood of a vertex is contained in the neighborhood of another vertex"""
        # different neighborhood
        for i, j in permutations(self.V, 2):
            # There must be a vertex adjecent to i which is not adjacent to j
            adjacentOnlyToI = []
            for k in self.V:
                if k == i or k == j:
                    continue
                kIsAdjacentOnlyToI = self.id()
                self.append([+self.var_edge(i, k), -kIsAdjacentOnlyToI])
                self.append([-self.var_edge(j, k), -kIsAdjacentOnlyToI])
                adjacentOnlyToI.append(kIsAdjacentOnlyToI)
            self.append([+self.var_edge(i, j)] + adjacentOnlyToI)

    def minConnectivity(self, connectivity_low) -> None:
        """Make graphs with minimum required connectivity

        :param connectivity_low: a lower bound on connectivity

        """
        assert self.n > connectivity_low  # an k-connected graph has at least k+1 vertices
        V = self.V
        var_edge = self.var_edge
        reachable = {
            (v, t, I): self.id() for k in range(args.connectivity_low) for I in combinations(sorted(set(V)), k) for v in set(V) - {min(set(V) - set(I))} - set(I) for t in V
        }  # u can reach v without I in t steps
        reachable_via = {
            (v, w, t, I): self.id()
            for k in range(args.connectivity_low)
            for I in combinations(sorted(set(V)), k)
            for v in set(V) - {min(set(V) - set(I))} - set(I)
            for t in V
            for w in set(V) - {min(set(V) - set(I)), v} - set(I)
        }  # u can reach v via w without I in t steps

        def var_reachable(v, t, I):
            return reachable[(v, t, I)]

        def var_reachable_via(v, w, t, I):
            return reachable_via[(v, w, t, I)]

        for k in range(args.connectivity_low):
            for I in combinations(sorted(set(V)), k):  # remove I and check if still connected
                u = min(set(V) - set(I))
                for v in set(V) - {u} - set(I):
                    for t in V:
                        if t == 0:
                            # reachable in first step if adjacent
                            self.append([-var_edge(v, u), +var_reachable(v, 0, I)])
                            self.append([+var_edge(v, u), -var_reachable(v, 0, I)])

                        else:
                            self.append([-var_reachable(v, t, I), +var_reachable(v, t - 1, I)] + [+var_reachable_via(v, w, t, I) for w in set(V) - set(I) - {v, u}])
                            self.append([+var_reachable(v, t, I), -var_reachable(v, t - 1, I)])
                            for w in set(V) - set(I) - {v, u}:
                                self.append([+var_reachable(v, t, I), -var_reachable_via(v, w, t, I)])

                                self.append([+var_reachable_via(v, w, t, I), -var_reachable(w, t - 1, I), -var_edge(w, v)])
                                self.append([-var_reachable_via(v, w, t, I), +var_reachable(w, t - 1, I)])
                                self.append([-var_reachable_via(v, w, t, I), +var_edge(w, v)])
                    # must be reached
                    self.append([+var_reachable(v, max(V), I)])

    def diameter2critical(self) -> None:
        """Ensure that the graph has diameter two and removing any edge results in a graph with diameter > 2"""
        g = self
        V = g.V
        var_edge = g.var_edge
        commonNeighbor = {(i, j, k): g.id() for i, j in combinations(V, 2) for k in set(V) - {i, j}}

        for i, j in combinations(V, 2):
            for k in set(V) - {i, j}:
                L = (i, j, k)
                g.CNF_AND_APPEND([+var_edge(i, k), +var_edge(j, k)], commonNeighbor[L])

        noCommonNeighbor = {(i, j): g.id() for i, j in combinations(V, 2)}
        for i, j in combinations(V, 2):
            for k in set(V) - {i, j}:
                # if the have a common neighbor, noCommonNeighbor is false
                g.append([-commonNeighbor[(i, j, k)], -noCommonNeighbor[(i, j)]])

        for i, j in combinations(V, 2):
            g.append([+var_edge(i, j)] + [+commonNeighbor[(i, j, k)] for k in set(V) - {i, j}])  # adjacent or common neighbor

        for i, j in combinations(V, 2):
            # ensure that critical i.e. if edge ij is present removing will lead to diamter > 2
            clause = [-var_edge(i, j), +noCommonNeighbor[(i, j)]]
            for k in set(V) - {i, j}:
                for v1, v2 in [(i, j), (j, i)]:
                    # v2 and k have diameter > after removing ij
                    # v1 adjacent to k and v1 is the only common neighbor from v2 and k. And k not adjacent to v2
                    diameterIncreasing = g.id()
                    g.append([+var_edge(v1, k), -diameterIncreasing])
                    g.append([-var_edge(v2, k), -diameterIncreasing])
                    for l in set(V) - {i, j, k}:
                        g.append([-commonNeighbor[(min(v2, k), max(v2, k), l)], -diameterIncreasing])
                    clause.append(diameterIncreasing)
            g.append(clause)

    def maxIndependentSet(self, x) -> None:
        """No independent sets of size greater than x

        :param x: upper bound on the size of independent sets

        """
        g = self
        for S in combinations(g.V, x + 1):
            g.append([+g.var_edge(i, j) for i, j in combinations(S, 2)])

    def maxClique(self, x) -> None:
        """No cliques of size greater than x

        :param x: upper bound on the size of cliques

        """
        g = self
        for S in combinations(g.V, x + 1):
            g.append([-g.var_edge(i, j) for i, j in combinations(S, 2)])

    def ckFree(self, k) -> None:
        """Forbid k-cycles (C_k) as non-induced subgraphs

        :param k: length of the cycle

        """
        g = self
        for cycle in permutations(g.V, k):
            if cycle[0] != min(cycle):
                continue
            if cycle[1] > cycle[-1]:
                continue
            g.append([-g.var_edge(cycle[i], cycle[(i + 1) % k]) for i in range(k)])  # at least one edge absent from potential cycle

    def minGirth(self, k) -> None:
        """Basic encoding to ensure that the girth is at least k, i.e., no cycle with length < k

        :param k: girth lower bound

        """
        for l in range(3, k):  # forbid 3 cycles up to k-1 cycles
            self.ckFree(l)

    def minGirthCompact(self, k) -> None:
        """Compact encoding to ensure that the girth is at least k, i.e., no cycle with length < k

        :param k: girth lower bound

        """
        g = self
        # check distance of i,j without edge i,j.
        for i, j in combinations(g.V, 2):
            reached = [g.var_edge(i, k) if k not in [i, j] and k > i else None for k in g.V]

            for _ in range(k - 4):  # if girth is 4 than no triangles so not in the loop
                reachedNew = [g.id() for _ in g.V]
                for k in g.V:
                    if k in [i, j] or k < i:
                        continue
                    g.append([-reached[k], +reachedNew[k]])  # already reached

                    # check if reached over l
                    for l in g.V:
                        if l in [i, j, k] or l < i:
                            continue
                        g.append([-g.var_edge(k, l), -reached[l], +reachedNew[k]])  # l reached in previous step and edge implies reached
                reached = reachedNew

            for k in g.V:
                if k in [i, j] or k < i:
                    continue
                # not reached, not adjacent to j, or edge not present
                g.append([-g.var_edge(i, j), -g.var_edge(j, k), -reached[k]])

    def maxChromaticNumber(self, chi) -> None:
        """The graph should be chi-colorable

        :param chi: nubmer of colors

        """
        g = self
        color = [[g.id() for _ in range(chi)] for _ in g.V]
        for v in g.V:
            g.append([color[v][i] for i in range(chi)])  # each vertex should have a color

        for v, w in combinations(g.V, 2):
            for i in range(chi):
                g.append([-g.var_edge(v, w), -color[v][i], -color[w][i]])  # adjacent vertices are not allowed to have the same color

    # --------------------------Encodings for initial static partition--------------------------------

    def sort_vertices_by_degree(self) -> None:
        """Create a static encoding of a partition, such that all vertices are sorted by degree.
        Currently only works for undirected graphs


        """
        g = self
        var_deg = []
        for v in g.V:
            counter_vars = g.counterFunction([g.var_edge(v, u) for u in g.V if u != v], g.n)
            var_deg.append(counter_vars)

        # compare consequent vertices
        for v1 in range(g.n - 1):
            v2 = v1 + 1

            # v2 is not allowed to have a lower degree
            for d in range(g.n - 1):
                g.append([-var_deg[v1][d], +var_deg[v2][d]])

            same_degree_options = []
            for d in range(g.n - 1):
                s = g.CNF_AND([+var_deg[v1][d], +var_deg[v2][d], -var_deg[v1][d + 1], -var_deg[v2][d + 1]])  # if s true then both have degree d - 1
                same_degree_options.append(s)
            g.CNF_OR_APPEND(same_degree_options, g.var_partition(v1, v2))


# ---------------------Main function------------------------------------------

if __name__ == "__main__":
    args, forwarding_args = getDefaultParser().parse_known_args()
    b = GraphEncodingBuilder(args.vertices, directed=args.directed, multiGraph=args.multigraph, staticInitialPartition=args.static_partition, underlyingGraph=args.underlying_graph, DEBUG=args.DEBUG)
    b.add_constraints_by_arguments(args)
    if args.no_solve:
        if args.cnf_file:
            with open(args.cnf_file, "w") as cnf_fh:
                b.print_dimacs(cnf_fh)
        else:
            b.print_dimacs(stdout)
    else:
        b.solveArgs(args, forwarding_args)
