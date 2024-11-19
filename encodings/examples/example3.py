from pysms.graph_builder import GraphEncodingBuilder
from itertools import combinations

b = GraphEncodingBuilder(7, directed=False)

chi = 3 # max chromatic number
color = [[b.id() for _ in range(chi)] for _ in b.V] # create ids; color[v][c] indicates whether v has color c
for v in b.V:
    b.append([color[v][i] for i in range(chi)])  # each vertex should have at least one color

for v, w in combinations(b.V, 2):
    for i in range(chi):
        b.append([-b.var_edge(v, w), -color[v][i], -color[w][i]])  # adjacent vertices are not allowed to have the same color

b.solve(allGraphs=True)
