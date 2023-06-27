# Applications



## Kochen-Specker vector systems and Kochen-Specker graphs

In the paper "Co-Certificate Learning with SAT Modulo Symmetries" we generated so called Kochen-Specker candidate graphs (short KS-candidates).
The script `./encodings/kochen_specker.py` can be used to reproduce the results. 

The following command can be used to generate all KS-candidates with 19 vertices:

```python
    python ./encodings/kochen_specker.py -v 19 --allGraphs
```

To generate cubes, we can use the arguments `--assignmentCutoffPrerunTime` which gives the time spend in seconds before cubing,  and `--assignmentCutoff`, which gives the minimal number of edge variables, which must be assigned for a cube.

```python
    python ./encodings/kochen_specker.py -v 22 --allGraphs --args_SMS "--assignmentCutoffPrerunTime 5 --assignmentCutoff 110"
```

The previous command produces cubes for 22 vertices.
To solve the problem with the cubes stored in the file `cubeFile`, we use the following command:

```python
    python ./encodings/kochen_specker.py -v 22 --allGraphs --args_SMS " --cubeFile cubeFile --cube2solve 2390 2392 "
```

## Erdős–Faber–Lovász conjecture

All encodings related to the EFL Conjecture are generated and solved by the script `./encodings/efl.py`.
The encoding is based on an incidence matrix to represent a hypergraph, i.e., the matrix indicates which vertices belong to which edge.

The most important arguments are as follows:

- `--n1 n`: set the number of vertices to n
- `--n2 m`: set the number of hyperedges to m
- `--chi3 i`: assert that the graph is (i-1)-reduced.
- `--maxClosedNeighborhood D`: ensure that the size of the neighborhood of a vertex is at most D and weakly covered.
- `--deactivateCovering`: deactivates the covering critierion; by default it is activated.

For example to test whether there is a hyper graph with 10 vertices and 14 edges, with edge chromatic number at least 11 which is 10-reduced,
we can use the following command

```python
    python ./encodings/efl.py --n1 10 --n2 14 --chi3 11 --args_SMS "--minEdgeChromaticNumberHypergraph 11"
```

Note that checking the minimum edge chromatic number is part of SMS and not the encoding itself and hence has to be given as argument to SMS.

Similar, we can check the FB Conjecture using the same script and the argument `maxClosedNeighborhood`.

```python
    python ./encodings/efl.py --n1 10 --n2 14 --maxClosedNeighborhood 7 --chi3 8 --args_SMS "--minEdgeChromaticNumberHypergraph 8"
```


## Planar graphs, Earth-Moon Problem, planar Turan Numbers and generation of OEIS integer sequences related to planar graphs

To enumerate planar graphs with different encodings, 
we used the script `./encodings/planar.py` with the following commands:

<ul>
  <li>Kuratowski based encoding:

```bash
python ./pysms/graph_builder.py -v $n --allGraphs --planar
```

    In this case, the planarity is not part of the encoding but rather forwarded to the SMS solver, and checked with a frequency of 1/5, i.e., only every 5th time, we check if the partially defined graph is planar. If not a suitable clause is added.
    </li>
  <li>
  Schnyder order based encoding:
```bash
python ./encodings/planarity.py -v $n --allGraphs --planar_schnyder
```
  </li>
  <li> Universal set based encoding:
```bash
python ./encodings/planarity.py -v $n --allGraphs --planar_universal
```
  </li>
</ul>



### Planar OEIS integer sequences

To investigate all OEIS sequences, we used the following commands:

$k$-connected $n$-vertex graphs 
for $k \in \{0,1,2,3,4,5\}$ (A88,A1349,A2218,A6290,A86216,A86217):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --connectivity_low $k
```

$k$-connected $n$-vertex planar graphs 
for $k \in \{0,1,2,3,4,5\}$ (A5470,A3094,A21103,A944,A7027,A361578):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --connectivity_low $k --planar
```

$k$-connected directed $n$-vertex graphs 
for $k \in \{0,1,2,3\}$ (A273,A3085,A361367,A361370):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --connectivity_low $k --directed
```

$k$-connected directed $n$-vertex planar graphs 
for $k \in \{0,1,2,3\}$ (A361366,A361368,A361369,A361371):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --connectivity_low $k --planar --directed
```

$k$-connected $n$-vertex triangulations
for $k \in \{3,4,5\}$ (A109,A7021,A111358):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --connectivity_low $k --planar --num_edges_low $((3*$n-6))
```

$n$-vertex planar graphs with even degrees (A49339):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --planar --evendegrees
```

connected $n$-vertex planar graphs with even degrees (A49365):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --planar --evendegrees --connectivity_low 1
```

$n$-vertex planar graphs with minimum degree at least $k$ 
for $k \in \{1,2,3,4,5\}$ (A49369-A49373):<br>
```bash
python ./pysms/graph_builder.py -v $n --allGraphs --planar --delta_low $k
```

connected triangle-free 3-regular $2n$-vertex planar graphs (A255600):<br>
```bash
python ./pysms/graph_builder.py -v $((2*$n)) --allGraphs --connectivity_low 1 --planar --Ckfree 3 --delta_low 3 --Delta_upp 3
```

2-connected 3-regular $n$-vertex planar graphs (A58378):<br>
```bash
python ./pysms/graph_builder.py -v $((2*$n)) --allGraphs --planar --delta_low 3 --Delta_upp 3 --connectivity_low 2
```


### Planar Turan numbers
For creating the encoding for the Turan numbers for planar graphs, we also use the script `./encodings/planarity.py`.
The arguments are as usual, especially we use the argument `--Ckfree c` to forbid all cycles with lentgth $c$.

For example, the following command produces a planar graph with $n$ vertices and atleast $m$ edges without a $4$-cycle or returns unsat.
```bash
python ./pysms/graph_builder.py -v $n --num_edges_low $m --Ckfree 4 --args_SMS " --planar 5 "
```


### Earth-Moon Problem


For creating the encoding for the Earth-Moon-Problem (based on planar directed graphs), we also use the script `./encodings/planarity.py` with the argument `--earthmoon c`, where $c$ is the minimum chromatic number of the searched for graph.

For example the following command produces a directed graph, whose underlying graph has $11$ vertices, chromatic number at least $9$ and is biplanar (i.e., has thickness $2$).
```bash
python ./encodings/planarity.py -v 11 --directed --earthmoon 9 --args_SMS " --thickness2 5"
```

Automatically, some assumptions are made when using the parameter `--earthmoon c`

- The graph $G_1$ is maximal planar
- $K_5$ and $K_{3,3}$ are excluded explizitly as subgraph for both $G_1$ and $G_2$.
- The underlying graph has minimum degree $\geq c - 1$.


Last, we can simply test whether the graph $C_5[4,4,4,4,3]$ is biplanar, using the following command
```bash
  python ./encodings/planarity.py -v 19 --directed --earthmoon_candidate1
```

