from pysms.graph_builder import GraphEncodingBuilder
from itertools import combinations

builder = GraphEncodingBuilder(7, directed=False)
for i, j, k in combinations(builder.V, 3):
    builder.append([-builder.var_edge(i, j), -builder.var_edge(i, k), -builder.var_edge(j, k)])
builder.solve(allGraphs=True)
