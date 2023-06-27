"""
    Script for building encodings related to incidence matrices.
    The underlying encoding is a bipartite graph.
"""

from pysms.graph_builder import *


def getParserIncidence():
    parser = getDefaultParser()
    # Remove required argument vertices
    for action in parser._actions:
        if "--vertices" in action.option_strings:
            parser._remove_action(action)  # remove vertex option
            break
    parser.add_argument("--n1", type=int, required=True, help="Number of rows of the incidence matrix")
    parser.add_argument("--n2", type=int, required=True, help="Number of columns of the incidence matrix")
    return parser


class IncidenceMatrixBuilder(GraphEncodingBuilder):
    """Incidence matrix represented as bipartite graph"""

    def __init__(self, n1, n2, staticInitialPartition=False):
        super().__init__(n=n1 + n2, staticInitialPartition=staticInitialPartition)
        del self.paramsSMS["vertices"]
        self.paramsSMS["bipartite"] = str(n1) + " " + str(n2)

        self.n1 = n1
        self.n2 = n2
        self.V1 = range(n1)
        self.V2 = range(n1, n1 + n2)

        # make bipartite
        for i, j in combinations(self.V1, 2):
            self.append([-self.var_edge(i, j)])
        for i, j in combinations(self.V2, 2):
            self.append([-self.var_edge(i, j)])

    def var_incident(self, i, j):
        return self.var_edge(i, j + self.n1)


if __name__ == "__main__":
    args = getParserIncidence().parse_args()
    b = IncidenceMatrixBuilder(args.n1, args.n2, staticInitialPartition=args.staticInitialPartition)
    b.add_constraints_by_arguments(args)
    b.solveArgs(args)
