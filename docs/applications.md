# Applications

## Kochen-Specker vector systems and Kochen-Specker graphs

In the paper [Co-Certificate Learning with SAT Modulo Symmetries](https://www.ijcai.org/proceedings/2023/0216.pdf) we generated so-called Kochen-Specker candidate graphs (short KS-candidates).
The script `./encodings/kochen_specker.py` can be used to reproduce the results. 

The following command can be used to generate all KS-candidates with 19 vertices:

```bash
    python ./encodings/kochen_specker.py -v 19 --all-graphs
```

To generate cubes, we can use the arguments `--assignment-cutoff-prerun-time` which gives the time spend in seconds before cubing,  and `--assignment-cutoff`, which gives the minimal number of edge variables, which must be assigned for a cube.

```bash
    python ./encodings/kochen_specker.py -v 22 --all-graphs --args-SMS "--assignment-cutoff-prerun-time 5 --assignment-cutoff 110"
```

The previous command produces cubes for 22 vertices.
To solve the problem with the cubes stored in the file `cubeFile`, we use the following command:

```bash
    python ./encodings/kochen_specker.py -v 22 --all-graphs --args-SMS " --cubes cubeFile --cube2solve 2390 2392 "
```

## Erdős–Faber–Lovász conjecture

All encodings related to the EFL Conjecture are generated and solved by the script `./encodings/efl.py`.
The encoding is based on an incidence matrix to represent a hypergraph, i.e., the matrix indicates which vertices belong to which edge. For all details we refer to the paper ([link](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SAT.2023.13)).

The most important arguments are as follows:

- `--n1 n`: set the number of vertices to n
- `--n2 m`: set the number of hyperedges to m
- `--selectCriticalSubgraph i`: assert that the graph is (i-1)-reduced.
- `--maxClosedNeighborhood D`: ensure that the size of the neighborhood of a vertex is at most D and weakly covered.
- `--deactivateCovering`: deactivates the covering critierion; by default it is activated.

For example to test whether there is a hyper graph with 10 vertices and 14 edges, with edge chromatic number at least 11 which is 10-reduced,
we can use the following command

```python
    python ./encodings/efl.py --n1 10 --n2 14 --selectCriticalSubgraph 11 --args-SMS "--min-chromatic-index-hypergraph 11"
```

Note that checking the minimum edge chromatic number is part of SMS and not the encoding itself and hence has to be given as argument to SMS.

Similar, we can check the FB Conjecture using the same script and the argument `maxClosedNeighborhood`.

```python
    python ./encodings/efl.py --n1 10 --n2 14 --deactivateCovering --maxClosedNeighborhood 7 --selectCriticalSubgraph 8 --args-SMS "--min-chromatic-index-hypergraph 8"
```


## Planar graphs, Earth-Moon Problem, planar Turan Numbers and generation of OEIS integer sequences related to planar graphs

To enumerate planar graphs, using different encodings, 
use `./encodings/planar.py` with the following commands (`$n` stands for the number of vertices). The theoretical part is described in the paper "SAT-Based Generation of Planar Graphs" ([link](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SAT.2023.14)):

<ul>
  <li>Kuratowski based encoding:

```bash
python -m pysms.graph_builder.py -v $n --all-graphs --planar
```

    In this case, the planarity is not part of the encoding but rather forwarded to the SMS solver, and checked with a frequency of 1/5, i.e., only every 5th time, we check if the partially defined graph is planar. If not a suitable clause is added.
    </li>
  <li>
  Schnyder order based encoding:
```bash
python ./encodings/planarity.py -v $n --all-graphs --planar_schnyder
```
  </li>
  <li> Universal set based encoding:
```bash
python ./encodings/planarity.py -v $n --all-graphs --planar_universal
```
  </li>
</ul>



### Planar OEIS integer sequences

To investigate all OEIS sequences, we used the following commands:

$k$-connected $n$-vertex graphs 
for $k \in \{0,1,2,3,4,5\}$ (A88,A1349,A2218,A6290,A86216,A86217):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --connectivity-low $k
```

$k$-connected $n$-vertex planar graphs 
for $k \in \{0,1,2,3,4,5\}$ (A5470,A3094,A21103,A944,A7027,A361578):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --connectivity-low $k --planar
```

$k$-connected directed $n$-vertex graphs 
for $k \in \{0,1,2,3\}$ (A273,A3085,A361367,A361370):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --connectivity-low $k --directed --underlying-graph
```

$k$-connected directed $n$-vertex planar graphs 
for $k \in \{0,1,2,3\}$ (A361366,A361368,A361369,A361371):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --connectivity-low $k --planar --directed --underlying-graph
```

$k$-connected $n$-vertex triangulations
for $k \in \{3,4,5\}$ (A109,A7021,A111358):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --connectivity-low $k --planar --num-edges-low $((3*$n-6))
```

$n$-vertex planar graphs with even degrees (A49339):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --planar --even-degrees
```

connected $n$-vertex planar graphs with even degrees (A49365):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --planar --even-degrees --connectivity-low 1
```

$n$-vertex planar graphs with minimum degree at least $k$ 
for $k \in \{1,2,3,4,5\}$ (A49369-A49373):<br>
```bash
python ./pysms/graph_builder.py -v $n --all-graphs --planar --delta-low $k
```

connected triangle-free 3-regular $2n$-vertex planar graphs (A255600):<br>
```bash
python ./pysms/graph_builder.py -v $((2*$n)) --all-graphs --connectivity-low 1 --planar --Ck-free 3 --delta-low 3 --Delta-upp 3
```

2-connected 3-regular $n$-vertex planar graphs (A58378):<br>
```bash
python ./pysms/graph_builder.py -v $((2*$n)) --all-graphs --planar --delta-low 3 --Delta-upp 3 --connectivity-low 2
```


### Planar Turan numbers
For creating the encoding for the Turan numbers for planar graphs, we also use the script `./encodings/planarity.py`.
The arguments are as usual, especially we use the argument `--Ck-free c` to forbid all cycles with lentgth $c$.

For example, the following command produces a planar graph with $n$ vertices and atleast $m$ edges without a $4$-cycle or returns unsat.
```bash
python ./pysms/graph_builder.py -v $n --num-edges-low $m --Ck-free 4 --args-SMS " --planar 5 "
```


### Earth-Moon Problem


For creating the encoding for the Earth-Moon-Problem (based on planar directed graphs), we also use the script `./encodings/planarity.py` with the argument `--earthmoon c`, where $c$ is the minimum chromatic number of the searched for graph.

For example the following command produces a directed graph, whose underlying graph has $11$ vertices, chromatic number at least $9$ and is biplanar (i.e., has thickness $2$).
```bash
python ./encodings/planarity.py -v 11 --directed --earthmoon 9 --args-SMS " --thickness2 5"
```

Automatically, some assumptions are made when using the parameter `--earthmoon c`

- The graph $G_1$ is maximal planar
- $K_5$ and $K_{3,3}$ are excluded explicitly as subgraph for both $G_1$ and $G_2$.
- The underlying graph has minimum degree $\geq c - 1$.


Last, we can simply test whether the graph $C_5[4,4,4,4,3]$ is biplanar, using the following command
```bash
  python ./encodings/planarity.py -v 19 --directed --earthmoon_candidate1
```


## Computing small Rainbow Cycle Numbers

For creating encodings related to computing the rainbow cycle number, we us the script `./encodings/efx.py`. The arguments `--partition-size` for the size of each block must be provided and 
also the number of vertices must be specified. The number of blocks is computed automatically by the number of vertices and the size of the blocks. 
For example 
```bash
python ./encodings/efx.py --directed --partition-size 3 -v 12
```
tries to compute a 4-partite graph without a rainbow cycle.
Adding the argument `--permutation` results in restricting the search to permutations and `--efx-propagator` ensure the absence of a rainbow cycle using a propagator instead of a static encoding.
For invariant pruning, we use a combination of two arguments for example `--delta-high-directed 6  --fix-first-vertex`. The first ensures that the outdegree is at most 6, the second that the first vertex cannot be permuted and it also must have outdegree exactly 6.


## QBF Problems ### {: #qbf}

In the AAAI'25 paper [Breaking Symmetries in Quantified Graph Search: A Comparative Study](), we evaluated several QBF (quantified Boolean formulas) solvers on a collection of quantified graph search problems: where the task is to generate graphs with coNP-hard properties such as non-colorability.
Below we explain how to install the solvers and generate the encodings from the paper.

### Solvers

The SMS package bundles a 2-QBF solver called 2Qiss.
It is available within the `smsg` binary.
Use `smsg --qcir-file <QCIR_FILE>` to pass a QCIR encoding, or call directly from `pysms.qcir_graph_builder`.

Several other solvers support SMS.
Follow the links and instructions found there to download and install them.
In order to use any of these solvers, SMS must first be installed.

- [Qfun](https://github.com/peitl/qfun-sms), pass the number of vertices with `-w`, `-E` to enumerate all graphs;
- [CQesto](https://github.com/peitl/cqesto-sms/tree/findcut) `-n` for the number of vertices, `-X` to enumerate. Note that CQesto enumerates over variables designated in a `free()` quantifier block;
- [Qute](https://github.com/peitl/qute-sms) `--sms-vertices` for the number of vertices, `--enumerate` to enumerate.

For each solver it is additionally to pass `--sms-cutoff` (equivalent to `smsg --cutoff`).

### Encodings

Below we list the exact PySMS commands used to generate the encodings.
The shell variables `VERTICES` and `QCIR_OUTPUT_FILE` are assumed to hold the number of vertices and the name of the file to write the encoding to, respectively.

### Coloring triangle-free graphs
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --chi-low 4 --Ck-free 3 --no-subsuming-neighborhoods --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --chi-low 5 --mtf --no-subsuming-neighborhoods --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
```

### Snarks
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --cubic --non-3-edge-colorable --Ck-free 3 --Ck-free 4 --two-connected --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
```

### Folkman graphs
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --folkmann 4 --no-subsuming-neighborhoods --print-qcir "$QCIR_OUTPUT_FILE"
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --folkmann 5 --no-subsuming-neighborhoods --print-qcir "$QCIR_OUTPUT_FILE"
```

### Kochen-Specker
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --kochen-specker --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
```

### Treewidth
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --tree-width-low-version2 4 --tree-width-upp-version1 4 --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --tree-width-low-version2 5 --tree-width-upp-version1 5 --all-graphs --print-qcir "$QCIR_OUTPUT_FILE"
```

### Domination in Cubic Graphs
```bash
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --domination-conjecture --connected-static --print-qcir "$QCIR_OUTPUT_FILE"
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --domination-conjecture --bipartite --print-qcir "$QCIR_OUTPUT_FILE"
python3 -m pysms.qcir_graph_builder -v "$VERTICES" --min-girth-compact 6 --domination-conjecture --print-qcir "$QCIR_OUTPUT_FILE"
```
