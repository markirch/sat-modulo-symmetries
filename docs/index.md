# SAT Modulo Symmetries

SAT Modulo Symmetries (SMS) is a framework that enhances SAT solvers with dynamic symmetry breaking. SMS is useful for isomorph-free generation and enumeration of graphs under constraints.

This SMS package contains the enhanced SAT solvers `smsg` and `smsd`, written in C++, which generate undirected and directed graphs modulo isomorphisms respectively, as well as a Python sub-module for easier creation of graph encodings with many predefined encodings to use.


## Requirements

- The SAT solver [Cadical](https://github.com/arminbiere/cadical).
- The library [Boost](https://www.boost.org/).
- To install the Python wrapper and encodings, [pip](https://pypi.org/project/pip).
- [CMake](https://cmake.org)

## Install

Currently, only installation on Linux is supported.

A copy of [CaDiCaL](https://github.com/arminbiere/cadical) is supplied (and expected) in the subdirectory `cadical`.
Build CaDiCaL with
```
cd cadical
./configure && make
```
as per its build instructions.

Next, install [Boost](https://www.boost.org/).
Then, SMS can be built and installed as follows: <!---CMake in the usual way, by executing the following commands: -->

```bash
git clone https://github.com/markirch/sat-modulo-symmetries
cd sat-modulo-symmetries
chmod +x build-and-install.sh
./build-and-install.sh
```

`build-and-install.sh` uses CMake, the built binaries are found in `build/src/`.
They are automatically installed to the default location (which probably requires root privileges, the script will ask).
Use `build-and-install.sh -l` instead if you don't have administrative priveleges and need to install locally (see `-h` for more options).

`smsg` generates undirected graphs, `smsd` is for directed graphs.

## Usage

SMS is a system that generates graphs modulo isomorphisms *under constraints*.
An example of such a task is to generate graphs with a certain maximum degree or number of edges, although such simple constraints barely scratch the surface of what SMS can do.

SMS consists of two main parts:

- the generation engine `smsg` (a compiled C++ binary) and its counterpart for directed graphs `smsd`; and
- a Python module, `pysms`, which comes with a library of pre-defined constraints, and an interface to implement custom new ones.

The typical workflow with SMS likewise has two steps.
The first is to generate the desired constraints and encode them in the DIMACS format.
This can be achieved using the module `pysms.graph_builder`, either selecting from the predefined constraints, or implementing a Python script to produce new ones.
The second is to 'solve' the constraints using `smsg` or `smsd`.
For convenience, the Python module has the capacity (turned on by default) to directly call `smsg` and solve the constraints it produces.
Therefore, for a quick start, we recommend using the Python wrapper.

### Python Module/Wrapper

The wrapper can be used as a library to build custom encodings or directly for enumerating graphs.
The wrapper is implemented in the file `pysms/graph_builder.py`, if you followed the installation process, the package `pysms` is installed and available for import. 
For example, the following command produces all graphs with 6 vertices up to isomorphism with minimum degree at least 3.

```bash
python -m pysms.graph_builder --vertices 6 --delta-low 3 --all-graphs
``` 

The found graphs are printed to standard output as Python lists of edges.

In order to only print the constraints (as a CNF formula in DIMACS) without solving, additionally pass `--no-solve`.

The most important options to `pysms.graph_builder` are as follows:

- `--vertices n` : search for graph with `n` vertices;
- `--all-graphs` : enumerate all graphs up to isomorphism satisfying the given properties (without this, the program terminates after finding the first graph);
- `--directed` : generate directed graphs (default is undirected);
- `--Delta_upp` : upper bound on the maximum degree;
- `--delta_low` : lower bound on the minimum degree;
- `--args_SMS s` : arguments forwarded to either `smsg` or `smsg`, i.e., the string `s` is appended as argument to the command line.

To get a complete list of all available options for the encoding builder, run 

```bash
python -m pysms.graph_builder --help
```

For all options available for `smsg` or `smsd`, run

```bash
smsg --help
```

The wrapper `pysms.graph_builder` will accept and forward any additional arguments to `smsg`.
There is the legacy option to pass arguments to `smsg` through `pysms.graph_builder --args-SMS="..."`, which is now deprecated.

SMS relies on a procedure called the _minimality check_ to filter out non-canonical isomorphic copies of graphs.
This procedure is often fast, but has worst-case exponential running time.
To avoid getting stuck in hard corner cases, we strongly recommend to use a time limit (cutoff) for the minimality check, by adding `--cutoff 20000`.
This limits the number of recursive calls in the minimality check, but potentially results in incomplete symmetry breaking.
The graphs can be filtered afterwards with tools like [Nauty](https://pallini.di.uniroma1.it/) using the `shortg` command.


#### Using the Encoding Builder in Python

The encoding builder is installed alongside `smsg` and `smsd` by `build-and-install.sh`, using `pip` as follows: 
```bash
pip install .
```
The installed importable module is `pysms`, but the pip package is called `sms-graph-builder`, and it can be removed with `pip uninstall sms-graph-builder` as usual.

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

The `builder` object contains the encoding and all necessary metadata such as the number of vertices and how to map edges to literals. 
The vertices are represented by the integers `0..n-1`.
A clause can be added using `builder.append(clause)`.  
The function `builder.var_edge(i,j)` returns the corresponding literal associated with the edge {i,j}.
So, `[-builder.var_edge(i,j), -builder.var_edge(i,k), -builder.var_edge(j,k)]` represents a clause that ensures that the vertex triple `i,j,k` doesn't form a triangle.

### Direct Usage

The most important arguments of `smsg` and `smsd` are as follows:

`smsg -v VERTICES [--all-graphs] [--frequency FREQUENCY] [--dimacs FILE] [--assignmentCutoff ACUTOFF] [--assignmentCutoffPrerunTime TIME] [--cutoff CUTOFF] [--printStats]`

- `-v VERTICES` gives the number of vertices (order) of the graph;
- `-b n m` searches for a bipartite graph with n+m vertices where the first n vertices and the last m vertices form the partitions;
- `--all-graphs` to compute all graphs;
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

- `--assignment-cutoff ACUTOFF` is used for cubing, if at least `ACUTOFF` edge variables are assigned and propagate is called, then a cube representing the partial assignment is generated.
- `--assignment-cutoff-prerun-time TIME` gives the prerun time in seconds before starting cubing.
- `--cutoff CUTOFF` gives the maximal recursive calls in the minimality check to avoid exponential behaviour. Note that this can result in an incomplete symmetry breaking.
- `--print-stats` prints some statistics of the computations, especially the time spend in different propagators.
- `--chi c` ensures chromatic number at least c.

For a complete list of all arguments call `smsg --help` and `smsd --help` respectively.

#### Examples

Generate all graphs on 7 vertices modulo isomorphism.

```bash
smsg -v 7 --all-graphs
```

Generate all non-bipartite graphs on 7 vertices (with chromatic number at least 3) and show times spend in minimality check and for coloring:

```bash
smsg -v 7 --all-graphs --chi 3 --print-stats
```

## Solvers

SMS currently requires the SAT solver [Cadical](https://github.com/arminbiere/cadical), but optionally also supports [Clingo](https://potassco.org/clingo). Please refer to the linked websites for installation instructions.
Cadical is pulled automatically as part of the build process.

### Installing Clingo

For the installation of Clingo see [https://potassco.org/clingo](https://potassco.org/clingo) and [https://github.com/potassco/clingo-cmake-example](https://github.com/potassco/clingo-cmake-example).
We recommend to install Clingo via conda, but other options should work as well.
