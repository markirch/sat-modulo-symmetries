from pysms.graph_builder import *

if __name__ == "__main__":
    parser = getDefaultParser()
    parser.add_argument("--partition-size", type=int, help="Size of the partition", required=True)
    parser.add_argument("--efx-propagator", action="store_true", help="Use the efx propagator")
    parser.add_argument("--efx-frequency", type=int, help="Frequency of the efx propagator", default=30)
    parser.add_argument("--additional-sym-breaking", action="store_true", help="Ensure that cannot map any element to larger element without introducing a cycle")
    parser.add_argument("--permutation", action="store_true", help="Ensure that mappings are permutations")
    parser.add_argument("--maximize-image", action="store_true", help="Try to maximize the image of the mappings without introducing additional cycles")
    parser.add_argument("--reachable-vertices-min", type=int, help="Minimum number of vertices, which can be reached by a colorful path (excluding the start vertex)")
    parser.add_argument("--fix-first-vertex", action="store_true", help="Fix the first vertex in the ordering and encodes that the reachable vertices are minimal")
    parser.add_argument("--vertex-ordering-file", type=str, help="File containing the vertex orderings for symmetry breaking", default="./vertexorderings.txt")
    parser.add_argument("--delta-high-directed", type=int, help="The maximal incoming degree")
    parser.add_argument("--delta-low-directed", type=int, help="The minimal incoming degree")

    parser.add_argument("--second-in-neighborhood-high", type=int, help="The maximal size of the second in-neighborhood (excluding vertices from the same block)")
    parser.add_argument("--second-out-neighborhood-high", type=int, help="The maximal size of the second in-neighborhood (excluding vertices from the same block)")
    parser.add_argument("--fix-first-second-neighborhood", action="store_true", help="Fix the first vertex for the previous two options if activated")

    # parser.add_argument("--pair-degree-high", type=int, help="The maximal incoming degree for pairs of vertices")
    parser.add_argument("--fix-first-pair", action="store_true")

    parser.add_argument("--max-indegree-from-other-partition", type=int, help="The maximal number of incomming edges from another partition")
    parser.add_argument("--fix-first-special", action="store_true", help="Fix the first vertex and the in-neighbors of the second partition")

    # parser.add_argument("--max-permutation-distance", type=int, help="For each equivalence class the maximal number of vertices not connected with a directed edge")

    args, forwarding_args = parser.parse_known_args()
    b = GraphEncodingBuilder(args.vertices, directed=args.directed, multiGraph=args.multigraph, staticInitialPartition=args.static_partition, underlyingGraph=args.underlying_graph, DEBUG=args.DEBUG)
    b.add_constraints_by_arguments(args)

    d = args.partition_size

    P = list(range(b.n // d))  # list of the partition labels
    partitioning = [list(range(i * d, (i + 1) * d)) for i in P]  # the concrete partitions

    # partition
    for p in partitioning:
        for i, j in permutations(p, 2):
            b.append([-b.var_edge_dir(i, j)])

    # for each vertex exactly one incomming edge from other blocks
    for p1, p2 in permutations(partitioning, 2):
        for v in p1:
            # at least one edge to others
            b.append([b.var_edge_dir(v, w) for w in p2])
            # at most one edge to others
            for w1, w2 in combinations(p2, 2):
                b.append([-b.var_edge_dir(v, w1), -b.var_edge_dir(v, w2)])

    # allRelevantDirectedEdgeVars = [b.var_edge_dir(i, j) for i, j in permutations(range(b.n), 2)]
    # testCount = 78
    # b.counterFunction(allRelevantDirectedEdgeVars, testCount, atLeast=testCount)

    if args.permutation:
        # ensure that permutation
        for p1, p2 in permutations(partitioning, 2):
            for v in p1:
                # at least one edge in other direction
                b.append([b.var_edge_dir(w, v) for w in p2])

    if False:
        # give satisfying assignment for testing
        for p1, p2 in permutations(P, 2):
            for v in range(d):
                if p1 < p2:
                    b.append([b.var_edge_dir(v + p1 * d, ((v + 1) % d) + p2 * d)])
                else:
                    b.append([b.var_edge_dir(v + p1 * d, (v % d) + p2 * d)])

    if args.efx_propagator:
        frequency = args.efx_frequency
        b.paramsSMS["efx"] = f"{len(partitioning)} {frequency} "
        b.paramsSMS["frequency"] = frequency
    else:

        # no colorful directed cycle, i.e., cycle such that each equivalence class is used only once
        colorFullPathVars = dict()
        for p1, p2 in permutations(P, 2):
            for v in partitioning[p1]:
                for w in partitioning[p2]:
                    otherPartitions = [p3 for p3 in P if p3 not in [p1, p2]]
                    for subset in chain.from_iterable(combinations(otherPartitions, r) for r in range(len(otherPartitions) + 1)):
                        colorFullPathVars[(v, w, frozenset(subset))] = b.id()

                    b.append([-colorFullPathVars[(v, w, frozenset(otherPartitions))], -b.var_edge_dir(w, v)])  # no color full path form v to w if there is a directed edge

        # ensure that colorFullPathVars are set to true if there is a path
        # - base case
        for p1, p2 in permutations(P, 2):
            for v in partitioning[p1]:
                for w in partitioning[p2]:
                    b.append([-b.var_edge_dir(v, w), +colorFullPathVars[(v, w, frozenset())]])  # there is a 1 path from v to w

        for p1, p2, p3 in permutations(P, 3):
            for v in partitioning[p1]:
                for w in partitioning[p2]:
                    for z in partitioning[p3]:
                        otherPartitions = [p4 for p4 in P if p4 not in [p1, p2, p3]]
                        for subset in chain.from_iterable(combinations(otherPartitions, r) for r in range(len(otherPartitions) + 1)):
                            b.append(
                                [-b.var_edge_dir(v, w), -colorFullPathVars[(w, z, frozenset(subset))], +colorFullPathVars[(v, z, frozenset(set(subset).union({p2})))]]
                            )  # if there is a path from v to w and a colorfol path from w to z than there is a colorful one from v to z but also using p2.

        # if there is a path then also for larger subsets (potentially better propagation)
        for p1, p2 in permutations(P, 2):
            for v in partitioning[p1]:
                for w in partitioning[p2]:
                    otherPartitions = [p3 for p3 in P if p3 not in [p1, p2]]
                    # print("otherPartitions", otherPartitions, file=sys.stderr)
                    for subset in chain.from_iterable(combinations(otherPartitions, r) for r in range(len(otherPartitions) + 1)):
                        for p3 in set(P) - {p1, p2} - set(subset):  # pick one partition that isn't in the subset and isn't p1 or p2
                            # when adding p3 to the subset and it there was a colorful path before, then there is still a colorful path
                            b.append(
                                [-colorFullPathVars[(v, w, frozenset(subset))], +colorFullPathVars[(v, w, frozenset(set(subset).union({p3})))]]
                            )  # if there is a path from v to w and a colorfol path from w to z than there is a colorful one from v to z but also using p2.

    def path(v, w, indicator_var):
        # There must be a directed path from v to w using each partition at most once if indicator_var is true
        for p1, p2 in permutations(P, 2):
            for v in partitioning[p1]:
                for w in partitioning[p2]:
                    pass
                    selectedVertices = [b.id() for _ in range(b.n)]  # indicating whether part of the path
                    b.append([selectedVertices[v]])
                    b.append([selectedVertices[w]])
                    # at most one vertex selected in each partition
                    for p in partitioning:
                        b.append([-selectedVertices[i], -selectedVertices[j]] for i, j in combinations(p, 2))

                    remainingEdges = [[None for _ in range(b.n)] for _ in range(b.n)]  # marking edges where both vertices are selected
                    for i, j in permutations(range(b.n), 2):
                        remainingEdges[i][j] = b.CNF_AND([selectedVertices[i], selectedVertices[j], b.var_edge_dir(i, j)])

                    # vhas exactly one out going edge
                    b.append([remainingEdges[v][j] for j in range(b.n) if j != v])
                    for j1, j2 in combinations(range(b.n), 2):
                        if v in [j1, j2]:
                            continue
                        b.append([-remainingEdges[v][j1], -remainingEdges[v][j2]])

                    # w has exactly one incomming edge
                    b.append([remainingEdges[j][w] for j in range(b.n) if j != w] + [-indicator_var])
                    for j1, j2 in combinations(range(b.n), 2):
                        if w in [j1, j2]:
                            continue
                        b.append([-remainingEdges[j1][w], -remainingEdges[j2][w]] + [-indicator_var])

                    # each other selected vertex has exactly one incomming and one outgoing edge
                    for i in range(b.n):
                        if i in partitioning[p1] or i in partitioning[p2]:
                            continue
                        b.append([remainingEdges[i][j] for j in range(b.n) if j != i] + [-indicator_var])
                        for j1, j2 in combinations(range(b.n), 2):
                            if i in [j1, j2]:
                                continue
                            b.append([-remainingEdges[i][j1], -remainingEdges[i][j2]] + [-indicator_var])

                        b.append([remainingEdges[j][i] for j in range(b.n) if j != i] + [-indicator_var])
                        for j1, j2 in combinations(range(b.n), 2):
                            if i in [j1, j2]:
                                continue
                            b.append([-remainingEdges[j1][i], -remainingEdges[j2][i]] + [-indicator_var])

    if args.reachable_vertices_min:
        useComplement = True if b.n - d - args.reachable_vertices_min < args.reachable_vertices_min else False  # use negation if cardinality is more compact to encode

        for p in P:
            for v in partitioning[p]:
                verticesInOtherPartition = [(w, p2) for p2 in P if p2 != p for w in partitioning[p2]]

                variablesIndicatingPath = []
                for w, p2 in verticesInOtherPartition:
                    otherPartitions = [p3 for p3 in P if p3 not in [p, p2]]
                    variablesIndicatingPath.append(colorFullPathVars[(v, w, frozenset(otherPartitions))])
                if useComplement:
                    variablesIndicatingPath = [-x for x in variablesIndicatingPath]
                    x = b.n - d - args.reachable_vertices_min
                    b.counterFunction(variablesIndicatingPath, x, atMost=x)
                else:
                    b.counterFunction(variablesIndicatingPath, args.reachable_vertices_min, atLeast=args.reachable_vertices_min)

        # for testing: assume first vertex to have the minimal neighbor of reachable vertices
        p = 0
        v = 0
        verticesInOtherPartition = [(w, p2) for p2 in P if p2 != p for w in partitioning[p2]]

        if args.fix_first_vertex:  # testing whether upperbounding speeds up the search
            variablesIndicatingPath = []
            for w, p2 in verticesInOtherPartition:
                otherPartitions = [p3 for p3 in P if p3 not in [p, p2]]
                variablesIndicatingPath.append(colorFullPathVars[(v, w, frozenset(otherPartitions))])
            if useComplement:
                variablesIndicatingPath = [-x for x in variablesIndicatingPath]
                x = b.n - d - args.reachable_vertices_min
                b.counterFunction(variablesIndicatingPath, x, atLeast=x)
            else:
                b.counterFunction(variablesIndicatingPath, args.reachable_vertices_min, atMost=args.reachable_vertices_min)

    if args.delta_high_directed:
        for p in partitioning:
            potentialNeighbors = [w for p2 in partitioning if p2 != p for w in p2]
            for v in p:
                b.counterFunction([b.var_edge_dir(w, v) for w in potentialNeighbors], args.delta_high_directed, atMost=args.delta_high_directed)

        if args.fix_first_vertex and not args.delta_low_directed:
            v = 0
            potentialNeighbors = [w for p2 in partitioning[1:] for w in p2]
            b.counterFunction([b.var_edge_dir(w, v) for w in potentialNeighbors], args.delta_high_directed, atLeast=args.delta_high_directed)

    if args.delta_low_directed:
        for p in partitioning:
            potentialNeighbors = [w for p2 in partitioning if p2 != p for w in p2]
            for v in p:
                b.counterFunction([b.var_edge_dir(w, v) for w in potentialNeighbors], args.delta_low_directed, atLeast=args.delta_low_directed)

        if args.fix_first_vertex:
            v = 0
            potentialNeighbors = [w for p2 in partitioning[1:] for w in p2]
            b.counterFunction([b.var_edge_dir(w, v) for w in potentialNeighbors], args.delta_low_directed, atMost=args.delta_low_directed)

    if args.max_indegree_from_other_partition:
        for p1, p2 in permutations(partitioning, 2):
            for v in p1:
                b.counterFunction([b.var_edge_dir(w, v) for w in p2], args.max_indegree_from_other_partition, atMost=args.max_indegree_from_other_partition)

        # fix first vertex to have maximum indegree to the second partition
        if args.fix_first_special:
            for j in range(args.max_indegree_from_other_partition):
                b.append([b.var_edge_dir(j + d, 0)])

    if args.second_in_neighborhood_high:
        for p in P:
            for v in partitioning[p]:
                # pass introduce variables which indicate for all other vertices whether they are in the second in-neighborhood
                secondNeighbors = []
                for p2 in P:
                    if p2 == p:
                        continue
                    for u in partitioning[p2]:
                        x = b.CNF_OR([ b.CNF_AND([b.var_edge_dir(u, w), b.var_edge_dir(w, v)]) for p3 in set(P) - {p, p2} for w in partitioning[p3]] + [b.var_edge_dir(u,v)]) # reachable over some other vertex w or directly
                        secondNeighbors.append(x)

                b.counterFunction(secondNeighbors, args.second_in_neighborhood_high, atMost=args.second_in_neighborhood_high)
                if v == 0:
                    secondNeighborsOfFirst = secondNeighbors

        if args.fix_first_second_neighborhood:
            b.counterFunction(secondNeighborsOfFirst, args.second_in_neighborhood_high, atLeast=args.second_in_neighborhood_high)

    if args.second_out_neighborhood_high:
        
        for p in P:
            for v in partitioning[p]:
                # pass introduce variables which indicate for all other vertices whether they are in the second in-neighborhood
                secondNeighbors = []
                for p2 in P:
                    if p2 == p:
                        continue
                    for u in partitioning[p2]:
                        x = b.CNF_OR([ b.CNF_AND([b.var_edge_dir(v, w), b.var_edge_dir(w, u)]) for p3 in set(P) - {p, p2} for w in partitioning[p3]] + [b.var_edge_dir(v,u)])  # reachable over some other vertex w
                        secondNeighbors.append(x)

                b.counterFunction(secondNeighbors, args.second_out_neighborhood_high, atMost=args.second_out_neighborhood_high)
                if v == 0:
                    secondNeighborsOfFirst = secondNeighbors
        
        if args.fix_first_second_neighborhood:
            b.counterFunction(secondNeighborsOfFirst, args.second_out_neighborhood_high, atLeast=args.second_out_neighborhood_high)

    partition = [d for _ in P]  # the partition sizes
    if True:
        allPermutations = True
        vertexOrderingsFile = args.vertex_ordering_file

        if allPermutations:
            # vertexOrderingsFile = f"./vertex_orderings{os.getpid()}_t{time.time()}.txt"
            with open(vertexOrderingsFile, "w") as f:
                
                if args.fix_first_special:
                    print(f"1 {partition[0] - 1} {args.max_indegree_from_other_partition} {d - args.max_indegree_from_other_partition} " + " ".join(map(str, partition[2:])), file=f)
                elif args.fix_first_vertex or args.fix_first_second_neighborhood:
                    # additionally fix first vertex
                    print(f"1 {partition[0] - 1} " + " ".join(map(str, partition[1:])), file=f)
                else:
                    print(" ".join(map(str, partition)), file=f)

                # swap pairs of partitions
                for permutedPartitioning in permutations(partitioning, len(partitioning)):
                    if args.fix_first_vertex and partitioning[0] != permutedPartitioning[0]:
                        continue

                    if  args.fix_first_special and (partitioning[0] != permutedPartitioning[0] or partitioning[1] != permutedPartitioning[1]):
                        continue # TODO

                    vertexOrdering = []
                    for p in permutedPartitioning:
                        vertexOrdering.extend(p)

                    print(" ".join(map(str, vertexOrdering)), file=f)
                    assert len(list(set(vertexOrdering))) == b.n
               
            b.paramsSMS["vertex-orderings-file"] = vertexOrderingsFile

        else:

            # vertexOrderingsFile = f"./vertex_orderings{os.getpid()}_t{time.time()}.txt"
            with open(vertexOrderingsFile, "w") as f:
                print(" ".join(map(str, partition)), file=f)

                print(" ".join(map(str, b.V)), file=f)

                # swap pairs of partitions
                for p1, p2 in combinations(P, 2):
                    vertexOrdering = []
                    for p in P:
                        if p == p1:
                            vertexOrdering.extend(partitioning[p2])
                        elif p == p2:
                            vertexOrdering.extend(partitioning[p1])
                        else:
                            vertexOrdering.extend(partitioning[p])

                    print(" ".join(map(str, vertexOrdering)), file=f)
                    assert len(list(set(vertexOrdering))) == b.n
            b.paramsSMS["vertex-orderings-file"] = vertexOrderingsFile

    else:
        b.paramsSMS["initial-partition"] = " ".join(map(str, partition))

    b.solveArgs(args, forwarding_args)
