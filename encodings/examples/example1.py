from pysms.graph_builder import GraphEncodingBuilder

builder = GraphEncodingBuilder(7, directed=False)
builder.maxDegree(3)
builder.solve(allGraphs=True)
