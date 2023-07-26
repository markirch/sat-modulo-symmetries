from pysms.incidince_matrix_builder import *

# we are looking for a linear hypergraph with
#   l  vertices
# whose intersection graph has
#   n  vertices
# and so the hypergraph has
#   n  hyperedges
#
# If the intersection graph is the primary object, we need
#   c(n, 2) variables for the intersection edges, and
#   l * n variables for the hypergraph incidence matrix
# and the intersection edges must come first for SMS.
#
# If the hypergraph is the primary object, we need
#   c(n+l, 2) variables for the hypergraph incidence matrix
#     including dummy vertex-vertex and hyperedge-hyperedge variables
#   c(n, 2) variables for the intersection edges
# and the hypergraph incidence variables must come first


class EFL(IncidenceMatrixBuilder):
    def __init__(self, n1, n2, staticInitialPartition=False):
        super().__init__(n1, n2, staticInitialPartition)

        self.vertices = n1
        self.m = n2

        E = range(self.m)
        V = range(self.vertices)
        N = range(self.m + self.vertices)  # all elements
        self.E = E
        self.VS = V  # self.V already used for underlying bipartite graph

        self.edge_vars_intersection = {(u, v): self.id() for u, v in combinations(E, 2)}
        self.paramsSMS["min-intersection-var"] = self.edge_vars_intersection[(0, 1)]  # set integer where intersection variables start
        self.paramsSMS["cutoff"] = 200000
        self.paramsSMS["hypergraph"] = ""
        self.paramsSMS["coloring-algo"] = 2  # use a CDCL solver for coloring

        self.hyperedges_share_vertex_v = [
            [[self.CNF_AND([self.edge_contains_vertex(e1, v), self.edge_contains_vertex(e2, v)]) if e1 < e2 else None for v in V] for e2 in E] for e1 in E
        ]  # hyperedges_share_vertex_v[e1][e2][v] true if v is in e1 and e2

        for e1, e2 in combinations(E, 2):
            # edge in intersection graph if one shared vertex
            self.CNF_OR_APPEND([+self.hyperedges_share_vertex_v[e1][e2][v] for v in V], self.var_edge_intersection_graph(e1, e2))
            for v1, v2 in combinations(V, 2):
                self.append([-self.hyperedges_share_vertex_v[e1][e2][v1], -self.hyperedges_share_vertex_v[e1][e2][v2]])  # ensures linearity, i.e., only one common vertex

        # minimum degree of a vertex 2
        for v in V:
            self.counterFunction([self.edge_contains_vertex(e, v) for e in E], 2, atLeast=2)

        # minimum size of a hyperedge 2
        for e2 in E:
            self.counterFunction([self.edge_contains_vertex(e2, v) for v in V], 2, atLeast=2)

        # covering
        # we can only assume this if we're not assuming vertex-criticality or different neighborhoods
        self.common_hyperedge = [[[self.id() if v1 < v2 else None for e in E] for v2 in V] for v1 in V]  # common_hyperedge[v1][v2][e] true if v1 and v2 are in e

        for v1, v2 in combinations(V, 2):
            for e2 in E:
                self.CNF_AND_APPEND([self.edge_contains_vertex(e2, v1), self.edge_contains_vertex(e2, v2)], self.common_hyperedge_func(v1, v2, e2))

    def common_hyperedge_func(self, v1, v2, e):
        return self.common_hyperedge[min(v1, v2)][max(v1, v2)][e]

    def edge_contains_vertex(self, e, v):
        return self.var_incident(v, e)  # !!!! different order; first vertices

    # variables to indicate whether an edge is present in the intersection graph
    def var_edge_intersection_graph(self, u, v):
        return self.edge_vars_intersection[(min(u, v), max(u, v))]

    def add_constraints_by_arguments(self, args):
        super().add_constraints_by_arguments(args)

        E = self.E
        V = self.VS

        # adjacent implies one with same label

        if not args.critical and not args.differentNeighborhood and not args.intersectionMinDegree and not args.deactivateCovering:
            for v1, v2 in combinations(V, 2):
                # some v has both incidence_graph
                self.append([self.common_hyperedge_func(v1, v2, v) for v in E])

        # Use hindman theorem stating that there must be at least 11 distinct vertices in hypereges with size > 2
        if args.hindman:
            # at least eleven incidence_graph which occur in a large edge
            selected = [self.id() for _ in V]
            selectedWitness = [[self.id() for _ in E] for _ in V]  # mark the large hyperedge with the label
            largeHyperEdge = []  # at least 3 incidence_graph
            for e2 in E:
                counter_vars = self.counterFunction([self.edge_contains_vertex(e2, v) for v in V], 3)
                largeHyperEdge.append(counter_vars[3 - 1])  # large edge has at least 3 vertices

            self.counterFunction(selected, 11, atLeast=11)
            for v in V:
                self.append([-selected[v]] + [selectedWitness[v][v] for v in E])
                for e2 in E:
                    self.append([-selectedWitness[v][e2], largeHyperEdge[e2]])
                    self.append([-selectedWitness[v][e2], self.edge_contains_vertex(e2, v)])

        # ensure an upperbound on the closed neighborhood of vertices in the hypergraph.
        if args.maxClosedNeighborhood:
            counter_vars_neighborhood = []

            adjacent_vertices = [[self.id() if v1 < v2 else None for v2 in V] for v1 in V]  # true if vertices occur in the same edge

            def adjacent_vertices_func(v1, v2):
                return adjacent_vertices[min(v1, v2)][max(v1, v2)]

            for v1, v2 in combinations(V, 2):
                clauses = CNF_OR([self.common_hyperedge_func(v1, v2, e) for e in E], adjacent_vertices_func(v1, v2))
                self.extend(clauses)

            if args.maxClosedNeighborhood > 1:
                # the normal case, limit the size of the neighborhood
                counter_range = args.maxClosedNeighborhood - 1
                upper_bound = counter_range
            else:
                # the combined case: allow any max closed neighborhood, but make sure edges cannot be added without raising it
                # and also request that the graph is colorable with fewer colors (this is to get non-extremal examples)
                counter_range = args.vertices - 1
                upper_bound = None

            for v in V:
                # counter over the open neighborhood
                counter_vars = self.counterFunction([adjacent_vertices_func(v, v2) for v2 in V if v2 != v], counter_range, atMost=upper_bound)

                if args.maxClosedNeighborhood > 1:
                    counter_vars_neighborhood.append(counter_vars[upper_bound - 1])  # check if contains exactly the maximum neighborhood
                else:
                    counter_vars_neighborhood.append(counter_vars)  # check if contains exactly the maximum neighborhood

            if args.maxClosedNeighborhood > 1:
                self.append(counter_vars_neighborhood)
                if args.neighborhoodMaximal:
                    for v1, v2 in combinations(V, 2):
                        self.append([adjacent_vertices[v1][v2], counter_vars_neighborhood[v1], counter_vars_neighborhood[v2]])
            else:
                # for each pair
                # the edge can be added if each endpoint has neighborhood size at most k-1
                # and there is another vertex that has neighborhood size at least k
                for v1, v2 in combinations(V, 2):
                    for k in range(1, args.vertices):
                        for v in V:
                            if v == v1 or v == v2:
                                continue
                            self.append([adjacent_vertices[v1][v2], counter_vars_neighborhood[v1][k - 1], counter_vars_neighborhood[v2][k - 1], -counter_vars_neighborhood[v][k - 1]])
                # self.append(counter_vars_neighborhood) # at least one vertex with highest possible neighborhood

                color = [[self.id() for c in V] for e in E]
                for e2 in E:
                    self.append(color[e2])
                for e1, e2 in combinations(E, 2):
                    for c in V:
                        self.append([-self.var_edge_intersection_graph(e1, e2), -color[e1][c], -color[e2][c]])

                for e2 in E:
                    self.append([-color[e2][args.vertices - 1]])
                    for c in range(1, args.vertices):
                        self.append([+counter_vars_neighborhood[v][c - 1] for v in V] + [-color[e2][c - 1]])

        # --------------------------------------arguments related to the chromatic index--------------------------------------
        if args.intersectionMinDegree:
            for i in E:
                self.counterFunction([self.var_edge_intersection_graph(i, j) for j in E if j != i], countUpto=args.intersectionMinDegree, atLeast=args.intersectionMinDegree)

        if args.intersectionMinDegreeLarge:
            for i in E:
                counter_intersections = self.counterFunction([self.var_edge_intersection_graph(i, j) for j in E if j != i], countUpto=args.intersectionMinDegreeLarge)
                counter_incidence_graph = self.counterFunction(
                    [self.edge_contains_vertex(i, v) for v in V],
                    countUpto=3,
                )
                # if at least 3 neighbors in the incidence graph (hyperedge size at least 3), then intersection degree at least as requested
                self.append([-counter_incidence_graph[3 - 1], +counter_intersections[args.intersectionMinDegreeLarge - 1]])

        if args.selectCriticalSubgraph:
            # select a subgraph which is critical
            selected = [self.id() for _ in E]
            self.counterFunction(selected, countUpto=args.selectCriticalSubgraph, atLeast=args.selectCriticalSubgraph)  # at least as many vertices as colors
            for i in E:
                counter_incidence_graph = self.counterFunction([self.edge_contains_vertex(i, v) for v in V], countUpto=3)
                # if at least 3 incidence_graph then degree then selected
                self.append([-counter_incidence_graph[3 - 1], +selected[i]])  # if >= 3 incidence_graph then definitely selected.

            for i in E:
                adjacent_and_selected = []
                for j in E:
                    if j == i:
                        continue
                    x = self.id()
                    clauses = CNF_AND([self.var_edge_intersection_graph(i, j), selected[j]], x)
                    self.extend(clauses)
                    adjacent_and_selected.append(x)
                counter_intersections = self.counterFunction(adjacent_and_selected, countUpto=args.selectCriticalSubgraph)
                # selected hyperedges must have at least selectCriticalSubgraph-1 many intersection neighbors
                self.append([-selected[i], +counter_intersections[args.selectCriticalSubgraph - 1 - 1]])

        # ------------------------------------not used currently------------------------------------------------------------------

        if args.vizing:
            maxDegree = args.vizing - 1
            for v in V:
                self.counterFunction([self.edge_contains_vertex(e, v) for e in E], countUpto=maxDegree, atMost=maxDegree)

        if args.minDegreeAtmost:
            # at least one vertex must have at most args.minDegreeAtmost
            selected = [self.id() for _ in E]
            self.append(selected)
            max_count = args.minDegreeAtmost + 1
            for e1 in E:
                counters = self.counterFunction([self.var_edge_intersection_graph(e1, e2) for e2 in E if e1 != e2], countUpto=max_count)
                self.append([-selected[e1], -counters[args.minDegreeAtmost + 1 - 1]])

                counters2 = self.counterFunction([self.edge_contains_vertex(e1, v) for v in V], countUpto=3)
                self.append([-selected[e1], counters2[3 - 1]])  # at least three vertices in that edge

        if args.differentNeighborhood:
            # different neighborhood
            for i, j in permutations(E, 2):
                # There must be a vertex adjecent to i which is not adjacent to j
                adjacentOnlyToI = []
                for k in E:
                    if k == i or k == j:
                        continue
                    kIsAdjacentOnlyToI = self.id()
                    self.append([+self.var_edge_intersection_graph(i, k), -kIsAdjacentOnlyToI])
                    self.append([-self.var_edge_intersection_graph(j, k), -kIsAdjacentOnlyToI])
                    adjacentOnlyToI.append(kIsAdjacentOnlyToI)
                self.append([+self.var_edge_intersection_graph(i, j)] + adjacentOnlyToI)

        if args.critical:
            # deletion of every hyperedge is colorable
            nColors = args.critical - 1
            for e2 in E:
                # check if G-v is args.critical - 1 colorable
                colors = [[self.id() for _ in E] for _ in range(nColors)]
                # at least one color
                for e1 in E:
                    if e1 != e2:
                        self.append([colors[r][e1] for r in range(nColors)])
                # adjacent once cannot have the same color
                for u1, u2 in combinations(E, 2):
                    if u1 == e2 or u2 == e2:
                        continue
                    for r in range(nColors):
                        self.append([-self.var_edge_intersection_graph(u1, u2), -colors[r][u1], -colors[r][u2]])

                # basic symmetry breaking for coloring
                # TODO smaller colors must not be available by previous vertices


if __name__ == "__main__":
    parser = getParserIncidence()
    parser.add_argument("--intersectionMinDegree", type=int, help="Request that each hyperedge intersects at least this many other hyperedges (holds whenever the intersection graph is critical).")
    parser.add_argument(
        "--intersectionMinDegreeLarge",
        type=int,
        help="Request that each large hyperedge (with at least 3 vertices) intersects at least this many other hyperedges (holds whenever the intersection graph contains a critical subgraph with all large hyperedges).",
    )
    parser.add_argument(
        "--selectCriticalSubgraph",
        type=int,
        help="Explicitly select a critical subgraph (for the given number of colors) of the intersection graph. The selected subgraph is required to contain all large hyperedges.",
    )
    parser.add_argument("--differentNeighborhood", help="Ensure that no hyperedge is a subset of another hyperedge", action="store_true")
    parser.add_argument("--critical", help="deleting any hyperedge decreases the chromatic index", action="store_true")
    # parser.add_argument('--primary', choices=["graph", "hypergraph"], help="choose which object SMS should operate on: the original linear hypergraph, or the intersection graph derived from it", default="graph")
    # parser.add_argument("--rowSwapStatic", action="store_true", help="Add a partial static symmetry breaking constraint that says that swapping any two consecutive rows (or columns) doesn't make the matrix lexicographically smaller")
    parser.add_argument("--hindman", action="store_true", help="The union of large hyperedges must have size at least 11 (Hindman proved in 1981 that otherwise the EFL conjecture holds)")
    parser.add_argument("--maxClosedNeighborhood", type=int, help="The maximum size of the closed neighborhood (= union of containing hyperedges) of a vertex in the hypergraph")
    parser.add_argument(
        "--neighborhoodMaximal",
        action="store_true",
        help="Request that the hypergraph be maximal with respect to closed neighborhood size (adding any hyperedge violates the neighborhood size constraint)",
    )
    parser.add_argument("--deactivateCovering", help="Admit non-covered hypergraphs (covered = each pair of vertices is contained in some hyperedge)", action="store_true", default=False)
    parser.add_argument("--vizing", type=int, help="Max degree + 1")
    parser.add_argument("--minDegreeAtmost", type=int, help="At least one vertex in the intersection graph must have degree at most this")
    args = parser.parse_args()
    b = EFL(args.n1, args.n2, staticInitialPartition=args.static_partition)
    b.add_constraints_by_arguments(args)
    b.solveArgs(args)
