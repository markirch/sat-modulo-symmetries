"""
Decodes models provided by SAT and QBF solvers in DIMACS notation
"""

from pysms.graph_builder import GraphEncodingBuilder

class Decoder:
    def __init__(self, n):
        self.n = n
        self.builder = GraphEncodingBuilder(n)

    def decode_line(self, line: str):
        lits = set(map(int, line.split()[1:]))
        G = []
        for u in range(self.n):
            for v in range(u+1, self.n):
                if self.builder.var_edge(u, v) in lits:
                    G.append((u, v))
        return G


if __name__ == "__main__":
    import sys
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-n", "--vertices", type=int, help="number of vertices")
    parser.add_argument("-x", "--extract", action="store_true", help="only output solutions, ignore other lines")
    args = parser.parse_args()

    dec = Decoder(args.vertices)
    for line in sys.stdin:
        if line.startswith("v"):
            print(dec.decode_line(line))
        elif not args.extract:
            print(line, end="")
