"""
    A script for checking the correctness of symmetry breaking clauses
"""

import argparse
from itertools import combinations

parser = argparse.ArgumentParser(description='Check the correctness of symmetry breaking clauses')
parser.add_argument('--vertices', type=int, help='The number of vertices in the graph', required=True)
parser.add_argument('--filename', type=str, help='The filename containing the clauses and witnessing permutation', required=True)
parser.add_argument('--directed', action='store_true', help='Whether the graph is directed or not', default=False)

args = parser.parse_args()


with open(args.filename, 'r') as f:
    for line in f:
        splittedLine = line.split(";")
        clause = splittedLine[0].strip()
        witness = splittedLine[1].strip()

        # clause is given as signed edge pairs for example (1,2) -(2,3) 
        # Split clause into individual edges, assuming spaces between them
        edges_with_signs = []
        edges = clause.split(" ")
        for edge in edges:
            if edge:  # Check if edge is not empty
                sign = '+' if edge[0] != '-' else '-'
                # Extract vertices assuming they are in the format (v1,v2) or -(v1,v2)
                vertices = edge.strip("()-").split(',')
                edge_tuple = (sign, (int(vertices[0]), int(vertices[1])))
                edges_with_signs.append(edge_tuple)

        # witness is given as a permutation of vertices
        witness = [int(vertex) for vertex in witness.split(" ")]

        # check whether witness is a permuation
        assert  sorted(witness) == list(range(args.vertices)) # must contain all vertices exactly once

        # check that the negation of the partially defined graph is indeed non-canonical 
        assert len(edges_with_signs) >= 2
        if not args.directed:
            for i,j in combinations(range(args.vertices), 2):
                
                pi = witness[i]
                pj = witness[j]
                if (i,j) == (pi, pj) or (i,j) == (pj, pi):
                    continue

                

                if len(edges_with_signs) == 2:
                    assert edges_with_signs[0][0] == "-"
                    assert edges_with_signs[0][1] == (i,j)

                    assert edges_with_signs[1][0] == "+"
                    assert edges_with_signs[1][1] == (pi, pj) or edges_with_signs[1][1] == (pj, pi)

                    break

                if edges_with_signs[0][0] == "-":
                    assert edges_with_signs[0][1] == (i,j) or edges_with_signs[0][1] == (j,i)
                else:
                    assert edges_with_signs[0][1] == (pi, pj) or edges_with_signs[0][1] == (pj, pi)

                edges_with_signs = edges_with_signs[1:]
        else:
            edgeOrdering = [] # The edge ordering used for the lexographical ordering
            for i in range(args.vertices):
                for j in range(i+1, args.vertices):
                    edgeOrdering.append((i,j))
                for j in range(i+1, args.vertices):
                    edgeOrdering.append((j,i))
            for i,j in edgeOrdering:
                pi = witness[i]
                pj = witness[j]
                if (i,j) == (pi, pj):
                    continue

                

                if len(edges_with_signs) == 2:
                    assert edges_with_signs[0][0] == "-"
                    assert edges_with_signs[0][1] == (i,j)


                    assert edges_with_signs[1][0] == "+"
                    assert edges_with_signs[1][1] == (pi, pj)
                    break

                if edges_with_signs[0][0] == "-":
                    assert edges_with_signs[0][1] == (i,j) 
                else:
                    assert edges_with_signs[0][1] == (pi, pj)

                edges_with_signs = edges_with_signs[1:]
                
            