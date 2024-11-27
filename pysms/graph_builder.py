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
import math


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
    solve_args.add_argument("--graph6-format", action="store_true", help="output graphs in graph6 format (Warning: Relatively slow)")
    # solve_args.add_argument("--graph-format", choices=["edge-list", "graph6"], default=None, help="output format of graphs. (Warning: Relatively slow)")

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
    constraint_args.add_argument("--color-critical", type=int, help="assert that deletion of any vertex makes the graph colorable with k-1 colors")

    constraint_args.add_argument("--fix-k-clique", type=int, help="the first k vertices are a clique (further it refines the initial partition)")

    constraint_args.add_argument("--circulant", action="store_true", help="ensure that graph is circulant and turnoff SMS")
    constraint_args.add_argument("--block-circulant", type=int, help="ensure that graph is block circulant with the given number of blocks and turnoff SMS")
    constraint_args.add_argument("--contains-cliques", type=int, nargs=2, help="ensures that the graph contains k cliques of size l and each edge is part of such a clique")
    constraint_args.add_argument("--partial-sym-break", action="store_true", help="break symmetry by considering swapping all pairs of vertices")
    constraint_args.add_argument("--bcp", action="store_true", help="apply bcp at the end of the encoding; not implemented very performant")
    constraint_args.add_argument("--exclude-balanced-bipartite-graph", action="store_true", help="exclude balanced bipartite graphs (only the once passing the static symmetry breaking)")

    constraint_args.add_argument("--fixed-induced-subgraph", type=str, help="file containing the fixed induced subgraph")
    constraint_args.add_argument("--fixed-induced-subgraph-line", type=int, help="choose line which is selected as induced subgraph")
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

        self.paramsSMS = {"vertices": self.n, "print-stats": True}  # default params

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

    def solve(self, allGraphs=False, hideGraphs=False, cnfFile=None, args_SMS="", forwarding_args=[], graph6_format=False) -> None:
        """Solve the formula, given the encoding, using SMS.

        :param allGraphs: Enumerate all satisfying graphs. Default value = False
        :param hideGraphs: Count all satisfying graphs. Default value = False
        :param cnfFile: Write constraints here, use a temporary if None. Default value = None
        :param forwarding_args: Forward these arguments to SMS. Default value = []
        :param graph_format: Toggle output format. Choices are [graph6, edge_list]. Default value = edge-list
        :param args_SMS: Forward these arguments to SMS. Deprecated, use `forwarding_args` instead. Default value = ""

        """
        if cnfFile == None:
            import time

            cnfFile = f"./temp{os.getpid()}_t{time.time()}.enc"  # TODO use tempfile module
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

        if graph6_format:
            import networkx as nx

            for line in os.popen(sms_command).read().split("\n"):
                if line and line[0] == "[":
                    edges = literal_eval(line)
                    print(nx.to_graph6_bytes(nx.Graph(edges), header=False).decode(), end="")
                elif self.DEBUG:
                    print(line, end="\n")
        else:
            os.system(sms_command)

        os.system(f"rm {cnfFile}")  # cleanup

    def solveArgs(self, args, forwarding_args) -> None:
        """Wrapper for solving using arguments provided by argsParser

        :param args: Arguments for PySMS
        :param forwarding_args: Arguments to be forwarded to SMS

        """
        self.solve(allGraphs=args.all_graphs, hideGraphs=args.hide_graphs, cnfFile=args.cnf_file, args_SMS=args.args_SMS, forwarding_args=forwarding_args, graph6_format=args.graph6_format)

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
        print(f"p cnf {self.nextId-1} {len(self)}", file=outstream)
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

        if args.color_critical:
            self.colorCritical(args.color_critical)

        if args.connectivity_low:
            self.minConnectivity(args.connectivity_low)

        if args.planar_kuratowski:
            self.paramsSMS["planar"] = 5  # DEFAULT planarity frequency

        if args.even_degrees:
            for u in self.V:
                shouldBe([+self.var_edge(u, v) for v in self.V if v != u], [i for i in self.V if i % 2 == 0], self, self, type=DEFAULT_COUNTER)

        if args.fix_k_clique:
            k = args.fix_k_clique
            for u, v in combinations(range(k), 2):
                self.append([self.var_edge(u, v)])
            self.paramsSMS["initial-partition"] = f"{k} {self.n - k}"

        if args.circulant:
            self.circulant()
        if args.block_circulant:
            self.block_circulant(args.block_circulant)

        if args.contains_cliques:
            self.contains_cliques(args.contains_cliques[0], args.contains_cliques[1])

        if args.partial_sym_break:
            self.partialSymBreak()

        if args.exclude_balanced_bipartite_graph:
            self.excludeBalancedBipartiteGraph()

        if args.fixed_induced_subgraph:
            with open(args.fixed_induced_subgraph, "r") as f:
                for i, line in enumerate(f):
                    if i == args.fixed_induced_subgraph_line:
                        l = line.split()
                        n = int(l[0])
                        l = l[1:]
                        edges = []
                        for i in range(len(l) // 2):
                            edges.append((int(l[2 * i]), int(l[2 * i + 1])))

                        matrix = [[0 for _ in range(n)] for _ in range(n)]
                        for u, v in edges:
                            matrix[u][v] = matrix[v][u] = 1
                        for u, v in combinations(range(n), 2):
                            if matrix[u][v] == 0:
                                self.append([-self.var_edge(u, v)])
                            else:
                                self.append([self.var_edge(u, v)])
                        self.paramsSMS["fixed-subgraph-size"] = n
                        break

        if args.bcp:  # !!!!!!!! must be applied at the end
            print("Number of propagated literals:", self.bcp(), file=stderr)

    # ------------degree encodings--------------

    def minDegree(self, delta, countertype=DEFAULT_COUNTER) -> list[list[int]]:
        """Minimum degree at least delta

        :param delta: degree lower bound
        :param countertype: specify a cardinality encoding. Default value = DEFAULT_COUNTER

        """
        return self.degreeBounds(self.V, delta, None, encoding=countertype)

    def maxDegree(self, delta, countertype=DEFAULT_COUNTER) -> list[list[int]]:
        """Maximum degree at most delta

        :param delta: degree upper bound
        :param countertype: specify a cardinality encoding. Default value = DEFAULT_COUNTER

        """
        return self.degreeBounds(self.V, None, delta, encoding=countertype)

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

    def minGirthCompact(self, minGirth) -> None:
        """
        Compact encoding to ensure that the girth is at least k, i.e., no cycle with length < k

        :param minGirth: girth lower bound
        """
        g = self
        # check distance of i,j without edge i,j.
        for i, j in combinations(g.V, 2):
            reached = [g.var_edge(i, k) if k not in [i, j] and k > i else None for k in g.V]

            for _ in range(minGirth - 4):  # if girth is 4 than no triangles so not in the loop
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

    def colorCritical(self, chi):
        nColors = chi - 1
        for v in self.V:
            # check if G-v is args.critical - 1 colorable
            colors = [[self.id() for _ in self.V] for _ in range(nColors)]
            # at least one color
            for u in self.V:
                if u != v:
                    self.append([colors[r][u] for r in range(nColors)])
            # adjacent once cannot have the same color
            for u1, u2 in combinations(self.V, 2):
                if u1 == v or u2 == v:
                    continue
                for r in range(nColors):
                    self.append([-self.var_edge(u1, u2), -colors[r][u1], -colors[r][u2]])

            # basic symmetry breaking: there must be a smaller vertex for each smaller color not choosen
            for v in self.V:
                for c in range(nColors):
                    for cSmaller in range(c):
                        self.append([-colors[c][v]] + [colors[cSmaller][u] for u in range(v)])

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

    def slimSingleStep(self, nOriginal, nReplaced, G, induced=False) -> None:
        """Apply single step of SLIM given the graph G as edgelist, where some random vertices are replaced"""
        import random

        VR = random.sample(range(nOriginal), k=nReplaced)

        nFixed = nOriginal - nReplaced
        self.paramsSMS["initial-partition"] = nFixed * "1 " + str(self.n - nFixed)

        remappingVertices = dict()
        cur = 0
        for v in set(range(nOriginal)) - set(VR):
            remappingVertices[v] = cur
            cur += 1
        r = remappingVertices

        remappedEdges = [(min(r[u], r[v]), max(r[u], r[v])) for u, v in G if u not in VR and v not in VR]
        # fix this subgraph in the encoding
        for u, v in remappedEdges:
            self.append([+self.var_edge(u, v)])
        if induced:
            for u, v in combinations(range(nFixed), 2):
                if (u, v) not in G and (v, u) not in G:
                    self.append([-self.var_edge(u, v)])

    def lex_smaller(self, seq1, seq2):
        """Ensure that seq1 is lexicographically smaller or equal than seq2"""
        assert len(seq1) == len(seq2)
        all_previous_equal = self.id()
        self.append([+all_previous_equal])
        for i in range(len(seq1)):
            self.append([-all_previous_equal, -seq1[i], +seq2[i]])  # all previous equal implies seq1[i] <= seq2[i]
            all_previous_equal_new = self.id()
            self.append([-all_previous_equal, -seq1[i], +all_previous_equal_new])
            self.append([-all_previous_equal, +seq2[i], +all_previous_equal_new])
            all_previous_equal = all_previous_equal_new

    def circulant(self, symmetryBreaking=True) -> None:
        self.paramsSMS["no-SMS"] = True
        for d in range(1, self.n // 2 + 1):  # math.ceil(self.n // 2) + 1): # math.ceil(self.n // 2) + 1
            for u in self.V:
                # must be equivalent
                u1 = 0
                u2 = d
                v1 = u
                v2 = (u + d) % self.n
                self.append([-self.var_edge(u1, u2), +self.var_edge(v1, v2)])
                self.append([+self.var_edge(u1, u2), -self.var_edge(v1, v2)])

        def edge2distance(u, v):
            return min(abs(u - v), self.n - abs(u - v))  # ensures that symmetric

        if False:
            # no-4-clique
            for potentialClique in combinations(range(1, self.n), 3):
                potentialClique = list(potentialClique) + [0]
                self.append([-self.var_edge(*e) for e in combinations(potentialClique, 2)])

        if False:
            deg = (self.n - 1) // 3 - 0
            self.counterFunction([self.var_edge(0, i) for i in range(1, self.n)], deg, atLeast=deg, counterType=DEFAULT_COUNTER)

        if False:
            # no-4-independent set
            for potentialClique in combinations(range(1, self.n), 3):
                potentialClique = list(potentialClique) + [0]
                self.append([+self.var_edge(*e) for e in combinations(potentialClique, 2)])

        # introducing any additional edge would result in an 4-clique
        if False:
            for u in range(1, self.n // 2 + 1):
                v = 0  # additional edge given by uv

                cliques = []
                for W in combinations(self.V, 2):
                    if u in W or v in W:
                        continue
                    clique = list(W) + [u, v]
                    cliquePairs = [p for p in combinations(clique, 2)]
                    distances = list(set(map(lambda p: edge2distance(p[0], p[1]), cliquePairs)))
                    distances.remove(edge2distance(v, u))  # remove the edges which would be introduced
                    # print(u, distances)
                    clique_var = self.CNF_AND([self.var_edge(0, d) for d in distances])
                    cliques.append(clique_var)
                # either edge is already present or it would introduce a 4-clique
                self.append([+self.var_edge(v, u)] + cliques)

        if symmetryBreaking:
            # symmetry breaking (only for primes)
            seq1 = [self.var_edge(0, i) for i in range(1, self.n // 2 + 1)]
            for first_vertex in [0]:  # don't need to check all because in same orbit
                for step_size in range(1, self.n // 2 + 1):
                    if first_vertex == 0 and step_size == 1:
                        continue

                    if math.gcd(step_size, self.n) != 1:
                        continue
                    reordering = [(first_vertex + step_size * i) % self.n for i in range(self.n)]
                    for i in range(self.n):
                        assert i in reordering
                    # print(reordering)

                    seq2 = [self.var_edge(reordering[0], reordering[i]) for i in range(1, self.n // 2 + 1)]
                    # seq2 = [ self.var_edge(0 , edge2distance(reordering[0], reordering[i])) for i in range(self.n // 2)]
                    self.lex_smaller(seq1, seq2)

    def block_circulant(self, nBlocks=None, symmetryBreaking=1) -> None:
        """if symmetry breaking = 0 then no sym breaking
        if symmetry brekaing = 1 then blocks
        if symmetry breaking = m then also permute m blocks
        # TODO also swap blocks
        """
        import sys

        print("Value of symmetry breaking", symmetryBreaking, file=sys.stderr)
        self.paramsSMS["no-SMS"] = True
        # Number of circulant graphs if not given then smallest prime factor
        if nBlocks == None:
            nBlocks = min([i for i in range(2, self.n) if self.n % i == 0])

        # ensure that the diagonal blocks are circulant
        blocksize = self.n // nBlocks
        for i in range(nBlocks):
            for j in range(i, nBlocks):
                for d in range(0 if i != j else 1, blocksize):  # also diagonal most be considered if i != j
                    for u in range(blocksize):
                        u1 = i * blocksize + 0
                        u2 = j * blocksize + d

                        v1 = i * blocksize + u
                        v2 = j * blocksize + (u + d) % blocksize

                        self.append([-self.var_edge(u1, u2), +self.var_edge(v1, v2)])
                        self.append([+self.var_edge(u1, u2), -self.var_edge(v1, v2)])

        if symmetryBreaking:
            # simple one only permuting each block and swap blocks; use vertex ordering such that blocks are sorted
            vertexPairOrdering = []
            print("Number of blocks:", nBlocks)
            for i in range(nBlocks):
                for u, v in combinations(range(i * blocksize, (i + 1) * blocksize), 2):
                    vertexPairOrdering.append((u, v))
            for u, v in combinations(self.V, 2):
                if (u, v) not in vertexPairOrdering:
                    vertexPairOrdering.append((u, v))
            symmetryCount = 0

            def checkGivenPermutation(p):
                seq1 = []
                seq2 = []
                for u, v in vertexPairOrdering:
                    up = p[u]
                    vp = p[v]
                    if {u, v} == {up, vp}:  # TODO eventuell consider duplicates
                        continue
                    seq1.append(self.var_edge(u, v))
                    seq2.append(self.var_edge(up, vp))
                self.lex_smaller(seq1, seq2)

            # swap blocks
            for i, j in combinations(range(nBlocks), 2):
                permutation = list(range(self.n))
                for k in range(blocksize):
                    permutation[blocksize * i + k] = blocksize * j + k
                    permutation[blocksize * j + k] = blocksize * i + k
                checkGivenPermutation(permutation)
                symmetryCount += 1

            print("Number of symmetries:", symmetryCount)

            # symmetry breaking on individual blocks
            for i in range(nBlocks):
                first_vertex = i * blocksize
                for step_size in range(2, blocksize // 2 + 1):
                    if math.gcd(step_size, blocksize) != 1:
                        continue
                    permutation = list(range(i * blocksize)) + [(step_size * i) % blocksize + first_vertex for i in range(blocksize)] + list(range((i + 1) * blocksize, self.n))
                    for i in range(self.n):
                        assert i in permutation
                    # print(reordering)
                    checkGivenPermutation(permutation)
                    symmetryCount += 1
            print("Number of symmetries:", symmetryCount)

        if symmetryBreaking and False:
            permutationsBlock = [[] for i in range(nBlocks)]  # for each block store the symmetries
            for i in range(nBlocks):
                first_vertex = i * blocksize
                for step_size in range(1, blocksize):
                    if math.gcd(step_size, blocksize) != 1:
                        continue
                    reordering = list(range(i * blocksize)) + [(step_size * i) % blocksize + first_vertex for i in range(blocksize)] + list(range((i + 1) * blocksize, self.n))
                    permutationsBlock[i].append(reordering)
            import sys

            # print(permutationsBlock, file=sys.stderr)
            counter = 0
            print("Value of symmetry breaking", symmetryBreaking, file=sys.stderr)
            for T in combinations(range(nBlocks), symmetryBreaking):
                print(T, file=sys.stderr)
                for P in product(*[permutationsBlock[i] for i in T]):
                    # print(P, file=sys.stderr)
                    counter += 1
                    seq1 = []
                    seq2 = []
                    for u, v in combinations(range(self.n), 2):
                        up = u
                        vp = v
                        for p in P:
                            up = p[up]
                            vp = p[vp]
                        if {u, v} == {up, vp}:
                            continue
                        seq1.append(self.var_edge(u, v))
                        seq2.append(self.var_edge(up, vp))
                    self.lex_smaller(seq1, seq2)
            print("Number of symmetries", counter, file=sys.stderr)
            # TODO also swap blocks

    def contains_cliques(self, s, k):
        """Ensures that the graph contains exactly s cliques of size k and each edge is part of such a clique"""
        cliqueVertices = [[self.id() for _ in self.V] for _ in range(s)]  # cliqueVertices[i][v] is true iff v is in clique i
        for i in range(s):
            self.counterFunction(cliqueVertices[i], k, k, k, counterType=DEFAULT_COUNTER)  # exactly k selected
            for u, v in combinations(self.V, 2):
                self.append([-cliqueVertices[i][u], -cliqueVertices[i][v], +self.var_edge(u, v)])
        # each edge is part of a clique
        for u, v in combinations(self.V, 2):
            self.append([-self.var_edge(u, v)] + [self.CNF_AND([cliqueVertices[i][u], cliqueVertices[i][v]]) for i in range(s)])

        # cliques can not be identical
        for i, j in combinations(range(s), 2):
            self.append([self.CNF_AND([+cliqueVertices[i][v], -cliqueVertices[j][v]]) for v in self.V])  # one in i but not in j

        # each vertex is in at most one clique
        for v in self.V:
            self.append([cliqueVertices[i][v] for i in range(s)])

        # symmetry breaking over cliques first clique is lex-smallest
        for i in range(1, s):
            self.lex_smaller(cliqueVertices[i - 1], cliqueVertices[i])

    def partialSymBreak(self):
        for u, v in combinations(self.V, 2):
            self.lex_smaller([self.var_edge(u, k) for k in range(self.n) if k not in [u, v]], [self.var_edge(v, k) for k in range(self.n) if k not in [u, v]])

    def bcp(self):
        """Perform unit propagation on the current encoding whilest ignoring the edge variables"""
        maxEdgeVar = self.n * (self.n - 1) // 2 if not self.directed else self.n**2 - self.n
        propagatedLiterals = [c[0] for c in self if len(c) == 1 if abs(c[0]) > maxEdgeVar]
        nPropagted = 0
        while len(propagatedLiterals) > 0:
            nPropagted += len(propagatedLiterals)
            # print("Propagated literals (inbetween):", len(propagatedLiterals))
            propagatedValues = [None for _ in range(self.nextId)]
            for l in propagatedLiterals:
                propagatedValues[abs(l)] = 1 if l > 0 else -1

            def sign(i):
                return 1 if i > 0 else -1

            def isSatisfied(c):
                return any([propagatedValues[abs(l)] and sign(l) == propagatedValues[abs(l)] for l in c])

            # already assumes that not satisfied
            def simplifyClause(c):
                return [l for l in c if not propagatedValues[abs(l)]]

            # delete literals from clauses and delete satisfied clauses
            newClauses = [simplifyClause(c) for c in self if not isSatisfied(c)]
            self.clear()
            self.extend(newClauses)
            propagatedLiterals = [c[0] for c in self if len(c) == 1 if abs(c[0]) > maxEdgeVar]

        return nPropagted

    def excludeBalancedBipartiteGraph(self):
        nHalf = [self.n // 2] if self.n % 2 == 0 else [self.n // 2, self.n // 2 + 1]
        for h in nHalf:
            self.append([self.var_edge(u, v) for u, v in combinations(range(h), 2)] + [self.var_edge(u, v) for u, v in combinations(range(h, self.n), 2)])


# ---------------------Main function------------------------------------------

if __name__ == "__main__":
    args, forwarding_args = getDefaultParser().parse_known_args()
    if forwarding_args:
        print("WARNING: Unknown arguments for python script which are forwarded to SMS:", forwarding_args, file=stderr)
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
