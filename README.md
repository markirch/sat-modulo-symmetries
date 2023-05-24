# SAT Modulo Symmetries

SAT Modulo Symmetries (SMS) is a framework that enhances SAT solvers with dynamic symmetry breaking. SMS is useful for isomorph-free generation and enumeration of graphs under constraints.

This SMS package contains the enhanced SAT solvers `smsg` and `smsd` which generate undirected and directed graphs modulo isomorphisms respectively, as well as a Python sub-module for easier creation of graph encodings.

## Solvers

SMS currently requires the SAT solver [Cadical](https://github.com/arminbiere/cadical), and optionally [Clingo](https://potassco.org/clingo). Please refer to the linked websites for installation instructions.
Cadical is pulled automatically as part of the build process.

SMS has the optional feature to generate graphs that do not contain any of a set of input forbidden subgraphs.
To enable, you need to install and build the [Glasgow Subgraph Solver](https://github.com/ciaranm/glasgow-subgraph-solver).
Run the script `pull-and-make-glasgow-solver.sh` to obtain and build the Glasgow Subgraph Solver, then configure with the flag `-DGLASGOW`.

## Installing Clingo

For the installation of Clingo see [https://potassco.org/clingo](https://potassco.org/clingo) and [https://github.com/potassco/clingo-cmake-example](https://github.com/potassco/clingo-cmake-example).
We recommend to install Clingo via conda, but other options should work as well.

## Install

With Cadical (and optionally Clingo) installed, SMS can be built with CMake. Execute the following commands:

```
git clone https://github.com/markirch/sat-modulo-symmetries
cd sat-modulo-symmetries
mkdir build
cmake -Bbuild -H.
cmake --build build
```

The built binaries will be found in `build/src/`. You can install them with

```
cmake --install build
```

`smsg` generates undirected graphs, `smsd` is for directed graphs.

To build with the Glasgow Subgraph Solver, configure with `cmake -Bbuild -H. -DGLASGOW=1` instead.

### Installing the Encoding Builder

To install the encoding builder run `pip install .`.
The installed pip package is called `sms-graph-builder`, and can be removed with `pip uninstall` as usual.

## Usage

The most important arguments of SMS are as follows:

`smsg -v VERTICES [--allModels] [--frequency FREQUENCY] [--cnf FILE] [--useCadical] [--assignmentCutoff ACUTOFF] [--assignmentCutoffPrerunTime TIME] [--cutoff CUTOFF] [--printStats]`

- `-v VERTICES` gives the number of vertices (order) of the graph;
- `-b n m` searches for a bipartite graph with $n+m$ vertices where the first $n$ vertices and the last $m$ vertices form the partitions;
- `--allModels` to compute all graphs;
- `--frequency FREQUENCY` used for balancing the time in the minimality check and the solver itself. For example if the frequency is 5 then the minimality check is only called every 5-th time (for reproducibility, note that this has a random component);
- `--dimacs FILE` the file name providing the SAT encoding in DIMACS format. For the undirected version the undirected edge variables are assumed to be given row-wise left-to-right from the upper triangle of the adjacency matrix. For the directed version the directed edge variables are assumed to be given row wise without the diagonal. For example, for an undirected graph on $n=4$ vertices the corresponding edge variables are given by the following literals: $e_{0,1} \leftrightarrow 1, e_{0,2} \leftrightarrow 2, e_{0,3} \leftrightarrow 3, e_{1,2} \leftrightarrow 4, e_{1,3} \leftrightarrow 5, e_{2,3} \leftrightarrow 6$. This can also be represented by the following matrix (ignoring diagonal entries): <br> 
`- 1 2 3 ` <br>
`1 - 4 5 ` <br>
`2 4 - 6 ` <br>
`3 5 6 - ` <br>
For the directed version the edge variables are in the following pattern: <br> 
`- 1 2 3 ` <br>
`4 - 5 6 ` <br>
`7 8 - 9 ` <br>
`10 11 12 - ` <br>
Note that this doesn't correspond with the order $\preceq$ of the vertex pairs used for defining the order on directed graphs.

- `--cnf FILE ` is the same as `--dimacs FILE` except assumes that each line (except lines starting with `c`) contains a clause without 0 at the end.
- `--useCadical` for using Cadical as solver. By default Clingo is used.
- `--assignmentCutoff ACUTOFF` is used for cubing, if at least `ACUTOFF` edge variables are assigned and propagate is called, then a cube representing the partial assignment is generated.
- `--assignmentCutoffPrerunTime TIME` gives the prerun time before starting cubing
- `--cutoff CUTOFF` gives the maximal recursive calls in the minimality check the avoid some bad cases. Note that this can result in an incomplete symmetry breaking.
- `--printStats` prints some statistics of the computations
- `-chi c` ensures chromatic number at least $c$.

## Examples

Generate all graphs on 7 vertices modulo isomorphism.

```smsg -v 7 --allModels```

Generate all non-bipartite graphs on 7 vertices (with chromatic number at least 3):

```smsg -v 7 --allModels -chi 3```

For more usage options see `src/main.cpp`

## Encoding Builder

To generate all graphs with 7 vertices and maximum degree at most 3, build and solve the encoding as follows:

```
from pysms.graph_builder import *
builder = GraphEncodingBuilder(7, directed=False)
builder.maxDegree(3)
builder.solve(allGraphs=True)
```
