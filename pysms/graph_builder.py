""" 
This file can either directly be used for generating all graphs described by some predefined properties which can be selected with flags,
or it can be use as starting point for constructing an encoding for graphs on your own.
"""

from itertools import *
from pysms.counters import *
import argparse
import os
from sys import *


def getDefaultParser():
    parser = argparse.ArgumentParser()  # default parser for some properties

    main_args = parser.add_argument_group(title="Main arguments", description="The number of vertices is mandatory, everything else is optional")
    solve_args = parser.add_argument_group(title="Solver options")
    constraint_args = parser.add_argument_group(title="Graph constraints", description="A set of pre-defined constraints for common applications, including applications from SMS papers")

    main_args.add_argument("--vertices", "-v", type=int, required=True, help="number of vertices")
    main_args.add_argument("--cnf-file", type=str, help="store the generated encoding here")
    main_args.add_argument("--directed", "-d", action="store_true", help="search for directed graphs")
    main_args.add_argument("--multigraph", "-m", type=int, help="search for a multigraph")
    main_args.add_argument("--underlying-graph", action="store_true", help="consider the underlying undirected graph for directed graphs")
    main_args.add_argument("--static-partition", action="store_true", help="specify a statically enforced partial vertex ordering (respected by SMS)")
    main_args.add_argument("--counter", choices=["sequential", "totalizer"], default="sequential", help="the CNF encoding for cardinality constraints")
    main_args.add_argument("--DEBUG", "-D", type=int, default=1, help="debug level")

    # solve options
    solve_args.add_argument("--all-graphs", "-a", action="store_true", help="generate all graphs (without this, the solver will exit after the first solution)")
    solve_args.add_argument("--hide-graphs", "-hg", action="store_true", help="do not display graphs (meant as a counting functionality, though the graphs still need to be enumerated)")
    solve_args.add_argument("--args-SMS", type=str, default="", help="command line to be appended to the call to smsg/smsd (see src/main.cpp or README.md)")

    # number of edges
    constraint_args.add_argument("--num-edges-upp", type=int, help="upper bound on the maximum number of edges")
    constraint_args.add_argument("--num-edges-low", type=int, help="lower bound on the minimum number of edges")

    # degree restrictions
    constraint_args.add_argument("--Delta-upp", type=int, help="upper bound on the maximum degree")
    constraint_args.add_argument("--delta-low", type=int, help="lower bound on the minimum degree")
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


def CNF_OR(ins, out):
    """Returns a list of clauses ensuring that `out` is equivalent to the disjunction of the literals in `ins`"""
    return [[-out] + ins] + [[out, -x] for x in ins]


def CNF_AND(ins, out):
    """Returns a list of clauses ensuring that `out` is equivalent to the conjunction of the literals in `ins`"""
    return [[out] + [-x for x in ins]] + [[-out, x] for x in ins]


DEFAULT_COUNTER = "sequential"


class GraphEncodingBuilder(IDPool, list):
    """A class for building an encoding in the context of SMS.
    The IDPool gives the next free variable whilst the list contains the clauses"""

    def __init__(self, n, directed=False, multiGraph=None, staticInitialPartition=False, underlyingGraph=False):
        super().__init__()
        self.directed = directed
        self.V = list(range(n))
        self.n = n
        self.DEBUG = 1
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
        """Get the propositional variable associated with the undirected edge {u,v}"""
        return self.varEdgeTable[u][v]

    def var_edge_multi(self, i, u, v) -> int:
        """Get the propositional variable associated with the undirected edge {u,v} of the i-th graph"""
        return self.varEdgeMultiTable[i][u][v]

    def var_edge_dir(self, u, v) -> int:
        """Get the propositional variable associated with the directed edge (u,v)"""
        return self.varEdgeDirectedTable[u][v]

    def var_partition(self, u, v) -> int:
        """Get the propositional variable which holds whether u and v are in the same partition"""
        return self.varStaticInitialPartition[u][v]

    def solve(self, allGraphs=False, hideGraphs=False, cnfFile=None, args_SMS=""):
        """Solve the formula, given the encoding, using SMS."""
        if cnfFile == None:
            cnfFile = f"./temp{os.getpid()}.enc"  # TODO use tempfile module
        self.print_dimacs(cnfFile)  # write script to temporary file

        program = "smsd" if self.directed else "smsg"  # we expect these binaries to be on PATH

        # add arguments

        if allGraphs:
            self.paramsSMS["all-graphs"] = ""
        if hideGraphs:
            self.paramsSMS["hide-graphs"] = ""
        if self.varStaticInitialPartition:
            self.paramsSMS["combine-static-dynamic"] = ""

        python_args_SMS = " ".join(f"--{param} {value}" for param, value in self.paramsSMS.items())

        sms_command = f"time {program} {python_args_SMS} {args_SMS} --dimacs {cnfFile}"  # TODO eventually parse args_SMS to allow to override

        if self.DEBUG:
            print("running the command: ", sms_command)
        stdout.flush()
        os.system(sms_command)
        os.system(f"rm {cnfFile}")

    def solveArgs(self, args):
        """Wrapper for solving using arguments provided by argsParser"""
        self.solve(allGraphs=args.all_graphs, hideGraphs=args.hide_graphs, cnfFile=args.cnf_file, args_SMS=args.args_SMS)

    # ------------------some utilies--------------------------

    def CNF_OR_APPEND(self, ins, out):
        self.extend(CNF_OR(ins, out))

    def CNF_AND_APPEND(self, ins, out):
        self.extend(CNF_AND(ins, out))

    def CNF_OR(self, ins):
        """returns a new variable which is true iff at least one of the literals in 'ins' is true"""
        out = self.id()
        self.extend(CNF_OR(ins, out))
        return out

    def CNF_AND(self, ins):
        """returns a new variable which is true iff all literals in 'ins' are true"""
        out = self.id()
        self.extend(CNF_AND(ins, out))
        return out

    def print_dimacs(self, filename=None):
        """Print the current encoding to the given file"""
        with open(filename, "w") as file:
            print(f"p cnf {len(self)} {self.nextId}", file=file)
            for c in self:
                print(" ".join(str(x) for x in c), 0, file=file)

    def counterFunction(self, variables, countUpto, atMost=None, atLeast=None, counterType="sequential"):
        """Wrapper for the counterFunction: constraints are added to the object itself and also the ids are given by the object."""
        return counterFunction(variables, countUpto, self, self, atMost=atMost, atLeast=atLeast, type=counterType)

    # -------------------------ENCODINGS----------------------------------------

    def add_constraints_by_arguments(self, args):
        """Add constraints based on args given by default parser"""
        g = self
        if g.DEBUG:
            print("Arguments:", args)
            stdout.flush()
        if args.Delta_upp:
            self.maxDegree(args.Delta_upp, args.counter)
        if args.delta_low:
            self.minDegree(args.delta_low, args.counter)
        if args.diam2_critical:
            self.diameter2critical()
        if args.num_edges_upp:
            self.numEdgesUpp(args.num_edges_upp, args.counter)
        if args.num_edges_low:
            self.numEdgesLow(args.num_edges_low, args.counter)

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

    def minDegree(self, delta, countertype=DEFAULT_COUNTER):
        """Minimum degree at least delta"""
        g = self
        for u in g.V:
            g.counterFunction([g.var_edge(u, v) for v in g.V if v != u], delta, atLeast=delta, counterType=countertype)

    def maxDegree(self, delta, countertype=DEFAULT_COUNTER):
        """Maximum degree at most delta"""
        g = self
        for u in g.V:
            g.counterFunction([g.var_edge(u, v) for v in g.V if v != u], delta, atMost=delta, counterType=countertype)

    # -------------------number of edges ------------------

    def numEdgesUpp(self, m, countertype=DEFAULT_COUNTER):
        """Upperbound on edges"""
        g = self
        g.counterFunction([g.var_edge(u, v) for u, v in combinations(g.V, 2)], m, atMost=m, counterType=countertype)

    def numEdgesLow(self, m, countertype=DEFAULT_COUNTER):
        """Lowerbound on edges"""
        g = self
        g.counterFunction([g.var_edge(u, v) for u, v in combinations(g.V, 2)], m, atLeast=m, counterType=countertype)

    # ------------------------------------------------------

    def mtf(self):
        self.ckFree(3)  # forbid triangles
        commonNeighbor = {(i, j, k): self.id() for i, j in combinations(self.V, 2) for k in self.V if k not in [i, j]}
        for i, j in combinations(self.V, 2):
            for k in self.V:
                if k in [i, j]:
                    continue
                self.CNF_AND_APPEND([self.var_edge(i, k), self.var_edge(j, k)], commonNeighbor[(i, j, k)])
            self.append([self.var_edge(i, j)] + [commonNeighbor[(i, j, k)] for k in self.V if k not in [i, j]])

    def noSubsumingNeighborhoods(self):
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

    def minConnectivity(self, connectivity_low):
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

    def diameter2critical(self):
        """Ensure that graph has diameter two and removing any edge results in a graph with diameter > 2"""
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

    def maxIndependentSet(self, x):
        """Maximal size of an independent set in the graph"""
        g = self
        for S in combinations(g.V, x + 1):
            g.append([+g.var_edge(i, j) for i, j in combinations(S, 2)])

    def maxClique(self, x):
        """Maximal size of a clique in the graph"""
        g = self
        for S in combinations(g.V, x + 1):
            g.append([-g.var_edge(i, j) for i, j in combinations(S, 2)])

    def ckFree(self, k):
        """Forbid k-cycles (C_k) as subgraphs"""
        g = self
        for cycle in permutations(g.V, k):
            if cycle[0] != min(cycle):
                continue
            if cycle[1] > cycle[-1]:
                continue
            g.append([-g.var_edge(cycle[i], cycle[(i + 1) % k]) for i in range(k)])  # at least one edge absent from potential cycle

    def minGirth(self, k):
        """Basic encoding to ensure that the girth is at least k, i.e., no cycle with length < k"""
        for l in range(3, k):  # forbid 3 cycles up to k-1 cycles
            self.ckFree(l)

    def minGirthCompact(self, k):
        """Compact encoding to ensure that the girth is at least k, i.e., no cycle with length < k"""
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

    def maxChromaticNumber(self, chi):
        """Maximum vertex chromatic number of the graph"""
        g = self
        color = [[g.id() for _ in range(chi)] for _ in g.V]
        for v in g.V:
            g.append([color[v][i] for i in range(chi)])  # each vertex should have a color

        for v, w in combinations(g.V, 2):
            for i in range(chi):
                g.append([-g.var_edge(v, w), -color[v][i], -color[w][i]])  # adjacent vertices are not allowed to have the same color

    # --------------------------Encodings for initial static partition--------------------------------

    def sort_vertices_by_degree(self) -> None:
        """
        Create by a static encoding a partition, such that all vertices are sorted.
        Currently only for undirected version
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
    args = getDefaultParser().parse_args()
    b = GraphEncodingBuilder(args.vertices, directed=args.directed, multiGraph=args.multigraph, staticInitialPartition=args.static_partition, underlyingGraph=args.underlying_graph)
    b.add_constraints_by_arguments(args)
    b.solveArgs(args)
