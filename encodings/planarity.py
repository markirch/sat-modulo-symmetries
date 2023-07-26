from pysms.graph_builder import *

"""
    This file contains some eager planarity encodings.
    Our experiments showed that kuratowski based planarity encodings are superiour,
    so we highly recommend using it.
"""


# planarity encodings
def planar_encoding_schnyder(V, var_edge, vpool, constraints, DEBUG=False):
    D = range(3)
    all_variables = []
    all_variables += [("u_smaller_v_i", (u, v, i)) for i in D for u, v in permutations(V, 2)]  # u < v in i-th linear u_smaller_v_i
    all_variables += [("uv_smaller_w_i", (u, v, w, i)) for i in D for u, v, w in permutations(V, 3)]  # u < w and v < w in i-th linear u_smaller_v_i

    all_variables_index = {}

    for v in all_variables:
        all_variables_index[v] = vpool.id()

    def var(L):
        return all_variables_index[L]

    def var_u_smaller_v_i(*L):
        return var(("u_smaller_v_i", L))

    def var_uv_smaller_w_i(*L):
        return var(("uv_smaller_w_i", L))

    if DEBUG:
        print("c\tdefine linear orders")
    for i in range(3):
        # anti-symmetrie + connexity
        for u, v in permutations(V, 2):
            constraints.append([+var_u_smaller_v_i(u, v, i), +var_u_smaller_v_i(v, u, i)])
            constraints.append([-var_u_smaller_v_i(u, v, i), -var_u_smaller_v_i(v, u, i)])

        # transitivity
        for u, v, w in permutations(V, 3):
            constraints.append([-var_u_smaller_v_i(u, v, i), -var_u_smaller_v_i(v, w, i), +var_u_smaller_v_i(u, w, i)])

    if DEBUG:
        print("c\tassert uv_smaller_w_i variable")
    for i in range(3):
        for u, v, w in permutations(V, 3):
            if u < v:
                constraints.append([-var_uv_smaller_w_i(u, v, w, i), +var_u_smaller_v_i(u, w, i)])
                constraints.append([-var_uv_smaller_w_i(u, v, w, i), +var_u_smaller_v_i(v, w, i)])
                constraints.append([+var_uv_smaller_w_i(u, v, w, i), -var_u_smaller_v_i(u, w, i), -var_u_smaller_v_i(v, w, i)])

    # if DEBUG: print("c\tantichain")
    for u, v in permutations(V, 2):
        constraints.append([+var_u_smaller_v_i(u, v, i) for i in D])

    if DEBUG:
        print("c\tplanarity criterion")
    # definition 1.1 from http://page.math.tu-berlin.de/~felsner/Paper/ppg-rev.pdf
    for u, v in combinations(V, 2):
        for w in set(V) - {u, v}:
            constraints.append([-var_edge(u, v)] + [var_uv_smaller_w_i(u, v, w, i) for i in range(3)])  # in upper triangle
            constraints.append([-var_edge(v, u)] + [var_uv_smaller_w_i(u, v, w, i) for i in range(3)])  # in lower triangle


def planarity_universal(V, var_edge, vpool, constraints):
    n = len(V)
    P_universal = {
        3: [(0, 0), (1, 0), (0, 1)],
        4: [(0, 0), (3, 0), (0, 3), (1, 1)],
        5: [(0, 3), (3, 3), (2, 2), (1, 2), (2, 0)],
        6: [(6, 0), (0, 0), (2, 2), (3, 2), (5, 3), (5, 6)],
        7: [(1, 10), (10, 2), (9, 2), (5, 3), (2, 8), (1, 9), (0, 0)],
        8: [(18, 0), (0, 0), (2, 2), (10, 7), (13, 5), (16, 6), (17, 6), (16, 18)],
        9: [(24, 15), (5, 0), (4, 5), (4, 14), (7, 15), (8, 16), (11, 16), (13, 19), (0, 24)],
        10: [(38, 0), (0, 92), (7, 90), (8, 88), (10, 85), (12, 79), (44, 8), (55, 22), (59, 27), (92, 66)],
        11: [(214, 0), (0, 13), (2, 16), (9, 26), (124, 12), (133, 11), (148, 9), (213, 1), (211, 4), (210, 6), (116, 179), (122, 197)],
    }

    # use minimal universal set or half grid if not available
    P = P_universal[n] if n in P_universal else [(x, y) for x in V for y in V if x < y]

    def sgn(x):
        return (x > 0) - (x < 0)

    def o3(a, b, c):
        (ax, ay) = a
        (bx, by) = b
        (cx, cy) = c
        return sgn((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))

    def pq_contains_r(p, q, r):
        dx1 = q[0] - p[0]
        dy1 = q[1] - p[1]
        dx2 = r[0] - p[0]
        dy2 = r[1] - p[1]
        assert dx1 != 0 or dy1 != 0
        if dx1 != 0 and (dx1 * dx2 < 0 or abs(dx2) > abs(dx1)):
            return False  # r not in segment pq
        if dy1 != 0 and (dy1 * dy2 < 0 or abs(dx2) > abs(dx1)):
            return False  # r not in segment pq
        return dx2 * dy1 == dx1 * dy2  # <=> r on line pq

    def pq_crosses_rs(p, q, r, s):
        return (o3(p, q, r) != o3(p, q, s)) and (o3(p, r, s) != o3(q, r, s)) or pq_contains_r(p, q, r) or pq_contains_r(p, q, s) or pq_contains_r(r, s, p) or pq_contains_r(r, s, q)

    PS = range(len(P))

    var_mapping = [[vpool.id() for p in P] for v in V]
    var_segment = [[vpool.id() if p1 < p2 else None for p2 in PS] for p1 in PS]
    # injective mapping
    for v in V:
        constraints.append([+var_mapping[v][p] for p in PS])
    for p in PS:
        for v1, v2 in combinations(V, 2):
            constraints.append([-var_mapping[v1][p], -var_mapping[v2][p]])
    for p1, p2 in combinations(PS, 2):
        for v in V:
            constraints.append([-var_mapping[v][p1], -var_mapping[v][p2]])

    # segments
    for u, v in combinations(V, 2):
        for p, q in permutations(PS, 2):
            constraints.append([-var_edge(u, v), -var_mapping[u][p], -var_mapping[v][q], +var_segment[min(p, q)][max(p, q)]])

    # no crossings
    for p, q, r, s in permutations(PS, 4):
        if pq_crosses_rs(P[p], P[q], P[r], P[s]):
            constraints.append([-var_segment[min(p, q)][max(p, q)], -var_segment[min(r, s)][max(r, s)]])


def getPlanarParser():
    parser = getDefaultParser()
    parser.add_argument("--planar_schnyder", "--ps", action="store_true", help="Planarity testing based on Schnyder orderings")
    parser.add_argument("--planar_universal", "--pu", action="store_true", help="Planarity testing based on universal point sets")
    parser.add_argument("--earthmoon", type=int, help="Criteria related to earth moon with given chromatic number; doesn't check thickness 2")
    parser.add_argument("--earthmoon_candidate1", action="store_true", help="Try to find decomposition of C5[4, 4, 4, 4, 3] into two planar graphs")
    parser.add_argument("--earthmoon_candidate2", action="store_true", help="Try to find decomposition of C7[4, 4, 4, 4, 4, 4, 4] into two planar graphs")
    return parser


class PlanarGraphBuilder(GraphEncodingBuilder):
    def __init__(self, n, directed=False, staticInitialPartition=False, underlyingGraph=False):
        super().__init__(n, directed, staticInitialPartition, underlyingGraph)

    def add_constraints_by_arguments(self, args):
        super().add_constraints_by_arguments(args)

        if args.planar_schnyder:
            planar_encoding_schnyder(self.V, self.var_edge, self, self)

        if args.planar_universal:
            planarity_universal(self.V, self.var_edge, self, self)

        if args.earthmoon:
            assert args.directed
            self.paramsSMS["frequency"] = 30
            V = self.V
            var_edge = self.var_edge_dir
            # The decomposition into two planar graphs is represented by double edges and edges in a single direction; w.l.o.g., one can assume that single edges are always at the lower triangular matrix.
            G = [[None for _ in V] for _ in V]
            for v, u in combinations(V, 2):
                G[v][u] = G[u][v] = var_edge(u, v)  # lower triangle gives union of both

            G1 = [[None for _ in V] for _ in V]
            for v, u in combinations(V, 2):
                G1[v][u] = G1[u][v] = var_edge(v, u)  # upper triangle matrix means at both

            G2 = [[None for _ in V] for _ in V]
            for v, u in combinations(V, 2):
                G2[v][u] = G2[u][v] = self.id()
                self.CNF_AND_APPEND([+var_edge(u, v), -var_edge(v, u)], G2[v][u])  # edge in G2 if single edge in lower half

            # explicitely ensure at least a certain minimum chromatic number
            chi = args.earthmoon
            from more_itertools import set_partitions

            for coloring in set_partitions(self.V, chi - 1):
                # print(x)
                clause = []  # monochromatic edge
                for color in coloring:
                    clause.extend([G[v][u] for v, u in combinations(color, 2)])
                self.append(clause)
            # self.paramsSMS["chi"] = args.earthmoon # set min chromatic number within propagator

            # only single edges in lower triangle matrix; i.e., if edge in upper than also in lower
            for v, u in combinations(V, 2):
                self.append([-var_edge(v, u), +var_edge(u, v)])  # if one upper half than also on lower half

            # w.l.o.g.: first graph is triangulation and second doesn't have more edges than triangulation
            if True:
                # euler and F=E*2/3 to get number of edges: V + F - E = 2; F = E * 2 / 3 -> V - 2= 1/3 E -> E = 3(V - 2)
                edgesInTriangulation = 3 * (len(V) - 2)
                self.counterFunction([G1[i][j] for i, j in combinations(V, 2)], countUpto=edgesInTriangulation, atMost=edgesInTriangulation, atLeast=edgesInTriangulation)
                self.counterFunction([G2[i][j] for i, j in combinations(V, 2)], countUpto=edgesInTriangulation, atMost=edgesInTriangulation)

            # graphs given by each direction (upper and lower triangular) must be planar
            if False:
                planar_encoding_schnyder(V, lambda v, u: G1[v][u], self, self)
                planar_encoding_schnyder(V, lambda v, u: G2[v][u], self, self)

            # mindegree on undirected version
            for i in V:
                self.counterFunction([G[i][j] for j in V if j != i], countUpto=args.earthmoon - 1, atLeast=args.earthmoon - 1)

            # forbid some trivial cases
            for A in combinations(V, 5):
                self.append([-G1[i][j] for i, j in combinations(A, 2)])
                self.append([-G2[i][j] for i, j in combinations(A, 2)])

            for A in combinations(V, 6):
                for B in combinations(A, 3):
                    if min(A) not in B:
                        continue
                    self.append([-G1[i][j] for i in set(A) - set(B) for j in B])
                    self.append([-G2[i][j] for i in set(A) - set(B) for j in B])

        if args.earthmoon_candidate1 or args.earthmoon_candidate2:
            partition = []
            E = []
            if args.earthmoon_candidate1:
                # C5[4, 4, 4, 4, 3]
                assert self.n == 19
                assert self.directed
                partition = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11], [12, 13, 14, 15], [16, 17, 18]]

            if args.earthmoon_candidate2:
                # C_7 where each vertex is replaced by K4s and joining neighborhoods
                assert self.n == 28
                assert self.directed
                for i in range(7):
                    partition.append(list(range(4 * i, 4 * (i + 1))))

            for P in partition:
                for i, j in combinations(P, 2):
                    E.append((i, j))
            for i in range(len(partition) - 1):
                for u in partition[i]:
                    for v in partition[i + 1]:
                        E.append((u, v))
            for u in partition[0]:
                for v in partition[-1]:
                    E.append((u, v))

            # print(E)
            for i, j in E:
                self.append([self.var_edge_dir(j, i)])

            for i, j in combinations(self.V, 2):
                if (i, j) not in E:
                    self.append([-self.var_edge_dir(i, j)])
                    self.append([-self.var_edge_dir(j, i)])

            self.paramsSMS["thickness2"] = "5"
            self.paramsSMS["initial-partition"] = " ".join(map(str, map(len, partition)))


args = getPlanarParser().parse_args()
b = PlanarGraphBuilder(args.vertices, directed=args.directed, staticInitialPartition=args.static_partition, underlyingGraph=args.underlying_graph)
b.add_constraints_by_arguments(args)
b.solveArgs(args)
