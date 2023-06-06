# SAT Modulo Symmetries

SAT Modulo Symmetries (SMS) is a framework that enhances SAT solvers with dynamic symmetry breaking. SMS is useful for isomorph-free generation and enumeration of graphs under constraints.

This SMS package contains the enhanced SAT solvers `smsg` and `smsd` written in C++ which generate undirected and directed graphs modulo isomorphisms respectively, as well as a Python sub-module for easier creation of graph encodings with many predefined encodings to use.



## Install

Currently, only installation on Linux is supported.
With [Cadical](https://github.com/arminbiere/cadical) installed, SMS can be built as follows: <!---CMake in the usual way, by executing the following commands: -->

```bash
git clone https://github.com/markirch/sat-modulo-symmetries
cd sat-modulo-symmetries
./buildScript.sh
```

<!--- Note that you can parallelize the build process, the precise way depends on your backend build system (for example with GNU Make do `cmake --build build -j n` to build with n concurrent processes). -->
The `buildScript.sh` uses CMake and the built binaries are found in `build/src/`. You can install them with

```bash
sudo cmake --install build
```

<!--- `cmake --install build` results in a permission denied on my laptop -->

`smsg` generates undirected graphs, `smsd` is for directed graphs.


## Usage

There are two main options to use the programs:

- Using a Python wrapper comming with many predefined properties which can be added 
- Using the programs `smsg`  and `smsd` directly.

For a quick start, we recommend using the Python wrapper.


### Python Wrapper

The wrapper can be used as a library to build custom encodings or directly used for enumerating graphs.
The wrapper is implemented in the file `pysms/graph_builder.py`. 
For example, the following command produces all graphs with 6 vertices up to isomorphism with minimum degree at least 3.

```bash
python pysms/graph_builder.py --vertices 6 --delta_low 3 --allGraphs
``` 
The graphs are printed to the output as edge list. The most relevant options are the following:

- `--vertices n` : search for graph with `n` vertices.
- `--allGraphs` : enumerate all graphs up to ismorphism, satisfying the given properties.
- `--directed` : generate directed graphs, per default undirected.
- `--Delta_upp` : upperbound on the maximum degree.
- `--delta_low` : lowerbound on the minimum degree.
- `--args_SMS s` : arguments forwarded to either `smsg` or `smsg`, i.e., string `s` is added as argument to the choosen program.

For a complete list of all available options for the python script, execute 
```bash
python pysms/graph_builder.py --help
```
For all options available to forward to `smsg` or `smsd`, execute
```bash
smsg --help
```


If too much time is spend in the minimality check, we strongly advice to use a cutoff of the minimality check algorithm, i.e., adding `--args_SMS "--cutoff 20000"`. This limits the number of recursive calls in the minimality check, but potentially results in an incomplete symmetry breaking. The graphs can be filtered afterwards with tools like  [nauty](https://pallini.di.uniroma1.it/) using the `shortg` command.


#### Installing the Encoding Builder

To install the encoding builder run 
```bash
pip install .
```
The installed pip package is called `pysms`, and can be removed with `pip uninstall` as usual.

#### Example 1

To generate all graphs with 7 vertices and maximum degree at most 3, build and solve the encoding as follows:

```python
from pysms.graph_builder import GraphEncodingBuilder
builder = GraphEncodingBuilder(7, directed=False)
builder.maxDegree(3)
builder.solve(allGraphs=True)
```

#### Example 2

The second example shows how one can enumerate all triangle free graphs with 7 vertices.

```python
from pysms.graph_builder import GraphEncodingBuilder
from itertools import combinations
builder = GraphEncodingBuilder(7, directed=False)
for i,j,k in combinations(builder.V, 3):
    builder.append([-builder.var_edge(i,j), -builder.var_edge(i,k), -builder.var_edge(j,k)])
builder.solve(allGraphs=True)
```

The `builder` object contains the encoding and all information like number of vertices and how to map edges to literals. 
A clause can be added using `builder.append(clause)`.  
The function `builder.var_edge(i,j)` returns the corresponding literal associated with the edge {i,j}.
So, `[-builder.var_edge(i,j), -builder.var_edge(i,k), -builder.var_edge(j,k)]` represents a clause that ensures that `i,j,k` doesn't form a triangle.

### Direct usage

The most important arguments of `smsg` and `smsd` are as follows:

`smsg -v VERTICES [--allGraphs] [--frequency FREQUENCY] [--dimacs FILE] [--assignmentCutoff ACUTOFF] [--assignmentCutoffPrerunTime TIME] [--cutoff CUTOFF] [--printStats]`

- `-v VERTICES` gives the number of vertices (order) of the graph;
- `-b n m` searches for a bipartite graph with n+m vertices where the first n vertices and the last m vertices form the partitions;
- `--allGraphs` to compute all graphs;
- `--frequency FREQUENCY` used for balancing the time in the minimality check and the solver itself. For example if the frequency is 5 then the minimality check is only called every 5-th time (for reproducibility, note that this has a random component);
- `--dimacs FILE` the file name providing the SAT encoding in DIMACS format. For the undirected version the undirected edge variables are assumed to be given row-wise left-to-right from the upper triangle of the adjacency matrix. For the directed version the directed edge variables are assumed to be given row wise without the diagonal. For example, for an undirected graph on n=4 vertices the corresponding edge variables are given by the following literals: e<sub>0,1</sub> &harr; 1, e<sub>0,2</sub> &harr; 2, e<sub>0,3</sub> &harr; 3, e<sub>1,2</sub> &harr; 4, e<sub>1,3</sub> &harr; 5, e<sub>2,3</sub> &harr; 6. This can also be represented by the following matrix (ignoring diagonal entries):

| | | | |
|---|---|---|---|
| - | 1 | 2 | 3 |
| 1 | - | 4 | 5 |
| 2 | 4 | - | 6 |
| 3 | 5 | 6 | - |

For the directed version the edge variables are given by the following pattern:

| | | | |
|---|---|---|---|
|  - |  1 |  2 |  3 |
|  4 |  - |  5 |  6 |
|  7 |  8 |  - |  9 |
| 10 | 11 | 12 |  - |

Note that this doesn't correspond with the order &preceq; of the vertex pairs used for defining the order on directed graphs.

- `--assignmentCutoff ACUTOFF` is used for cubing, if at least `ACUTOFF` edge variables are assigned and propagate is called, then a cube representing the partial assignment is generated.
- `--assignmentCutoffPrerunTime TIME` gives the prerun time in seconds before starting cubing.
- `--cutoff CUTOFF` gives the maximal recursive calls in the minimality check to avoid exponential behaviour. Note that this can result in an incomplete symmetry breaking.
- `--printStats` prints some statistics of the computations, especially the time spend in different propagators.
- `--chi c` ensures chromatic number at least c.

For a complete list of all arguments call `smsg --help` and `smsd --help` respectively.

#### Examples

Generate all graphs on 7 vertices modulo isomorphism.

```bash
smsg -v 7 --allGraphs
```

Generate all non-bipartite graphs on 7 vertices (with chromatic number at least 3) and show times spend in minimality check and for coloring:

```bash
smsg -v 7 --allGraphs --chi 3 --printStats
```

## Solvers

SMS currently requires the SAT solver [Cadical](https://github.com/arminbiere/cadical), but optionally als supports [Clingo](https://potassco.org/clingo). Please refer to the linked websites for installation instructions.
Cadical is pulled automatically as part of the build process.

### Installing Clingo

For the installation of Clingo see [https://potassco.org/clingo](https://potassco.org/clingo) and [https://github.com/potassco/clingo-cmake-example](https://github.com/potassco/clingo-cmake-example).
We recommend to install Clingo via conda, but other options should work as well.