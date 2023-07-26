from pysms.graph_builder import *
from itertools import *


class KSGraphEncodingBuilder(GraphEncodingBuilder):
    def __init__(self, n, staticInitialPartition=False):
        super().__init__(n, staticInitialPartition=staticInitialPartition)

        # Create triangle variables needed for encoding
        self.triangleVariables = {(i, j, k): self.id() for i, j, k in combinations(self.V, 3)}
        for i, j, k in combinations(self.V, 3):
            self.CNF_AND_APPEND([self.var_edge(i, j), self.var_edge(j, k), self.var_edge(i, k)], self.triangleVariables[(i, j, k)])  # ensure that true iff it is a triangle

        self.paramsSMS["triangle-vars"] = self.triangleVariables[(0, 1, 2)]
        self.paramsSMS["non010"] = ""
        self.paramsSMS["frequency"] = 30
        # self.paramsSMS["forbiddenSubgraphs"] = "./unembeddableSmall.txt"

        # ensure properties of KS graph
        self.maxChromaticNumber(4)  # at most chromatic number 4
        self.minDegree(3)  # minimum degree 3
        self.ckFree(4)  # no 4 cycle
        self.eachVertexInTriangle()

    def eachVertexInTriangle(self):
        g = self
        for v in g.V:
            g.append([g.triangleVariables[(i, j, k)] for i, j, k in combinations(g.V, 3) if v in [i, j, k]])  # check that at least one triangle is present

    def eachEdgeInTriangle(self):
        g = self
        for u,v in combinations(g.V,2):
            g.append([-g.var_edge(u,v)] + [g.triangleVariables[(i, j, k)] for i, j, k in combinations(g.V, 3) if v in [i, j, k] and u in [i, j, k]])  # check that at least one triangle is present

    def add_constraints_by_arguments(self, args):
        super().add_constraints_by_arguments(args)
        if args.edge_in_triangle:
            self.eachEdgeInTriangle()


if __name__ == "__main__":
    parser = getDefaultParser()
    parser.add_argument("--edge_in_triangle", action="store_true", help="Ensure that each edge is in a triangle")
    args = parser.parse_args()
    g = KSGraphEncodingBuilder(args.vertices, staticInitialPartition=args.static_partition)
    g.add_constraints_by_arguments(args)
    g.solveArgs(args)
