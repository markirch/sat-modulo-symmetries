# SAT Modulo Symmetries

SAT Modulo Symmetries (SMS) is a framework and software package for isomorph-free generation and enumeration of graphs under constraints.
Constraints for SMS can be specified declaratively, in a combination of propositional logic ([DIMACS](https://jix.github.io/varisat/manual/0.2.0/formats/dimacs.html)) and [custom extensions](advanced), and solved using specially adapted [SAT solvers](https://en.wikipedia.org/wiki/SAT_solver).
With SMS, you can easily implement and solve complex graph constraints thanks to the combination of a flexible input language and a powerful SAT solver.


SMS contains:

- the enhanced solvers `smsg` and `smsd`, based on the award-winning SAT solver [CaDiCaL](https://github.com/arminbiere/cadical), which are the engine that generates undirected and directed graphs modulo isomorphisms;
- C++ SMS libraries (`libsms_static.a`, `libsmsdir_static.a`, `libsms.so`, `libsmsdir.so`);
- the Python library `PySMS`, which contains a collection of pre-defined constraints and an interface to start building your own graph generation scripts;
- a collection of examples and applications;
- a 2QBF solver for modeling graphs with coNP properties. 

## Requirements

In order to build SMS, you will need the following:

- a C++20-ready compiler;
- [Boost](https://www.boost.org/) 1.74 or newer, including `libboost-program-options-dev`;
- [CMake](https://cmake.org) 3.12 or newer; and
- to install PySMS, [pip](https://pypi.org/project/pip).

A copy of [CaDiCaL](https://github.com/arminbiere/cadical) is supplied with SMS.

## Install

Currently, only installation from source on Linux is officially supported.

Make sure you have [Boost](https://www.boost.org/) installed.
A copy of [CaDiCaL](https://github.com/arminbiere/cadical) is supplied (and expected) in the subdirectory `cadical`.
Clone and enter the repository

```bash
git clone https://github.com/markirch/sat-modulo-symmetries
cd sat-modulo-symmetries
```

and build SMS with:

```bash
chmod +x build-and-install.sh
./build-and-install.sh
```

`build-and-install.sh` uses CMake, the built binaries are found in `build/src/`.
The executables, static libraries, and required headers are automatically installed to the default location (something like `/usr/local/`, which probably requires root privileges, for which `build-and-install.sh` will ask).

!!! note
	Use `build-and-install.sh -l` instead if you don't have administrative priveleges and need to install locally (see `-h` for more options).

!!! warning
	An earlier version of this manual recommended to first build CaDiCaL. This step is now included in `build-and-install.sh` and **should not** be performed separately. If you build CaDiCaL separately, you need to make sure it is configured with `-fPIC`.

`smsg` generates undirected graphs, `smsd` is for directed graphs.
Likewise use `libsms_static.a`, `libsms.so` to generate undirected graphs, and `libsmsdir_static.a`, `libsmsdir.so` for directed graphs.

## Troubleshooting

If the installation of PySMS fails, try updating pip with `python3 -m pip install --update pip`.

You may see the `Illegal instruction` error when running SMS built with the Glasgow Subgraph Solver (see [advanced usage](http://localhost:8000/advanced/#forbidden-subgraphs)) on a different machine than where it was built.
If this happens, comment out (with `#`) the lines in `glasgow-subgraph-solver/CMakeLists.txt` that talk about `-march=native`, delete `glasgow-subgraph-solver/build`, and re-run `./build-and-install.sh -s`.

Please report any other issues [to us by email](team) or open an issue on [GitHub](https://github.com/markirch/sat-modulo-symmetries/issues).

## Overview

SMS is a system that generates graphs modulo isomorphisms *under constraints*.
An example of such a task is to generate graphs with a certain maximum degree or number of edges, although such simple constraints barely scratch the surface of what SMS can do.

SMS consists of two main parts:

- the generation engine `smsg` (a compiled C++ binary) and its counterpart for directed graphs `smsd` (and the corresponding static and shared libraries); and
- a Python module, PySMS, which comes with a library of pre-defined constraints, and an interface to implement custom new ones.

The typical workflow with SMS likewise has two steps.
The first is to generate the desired constraints and encode them in the [DIMACS](https://www.cs.utexas.edu/users/moore/acl2/manuals/current/manual/index-seo.php/SATLINK____DIMACS) format.
This can be achieved using the module `pysms.graph_builder`, either selecting from the predefined constraints, or implementing custom ones.
The second is to 'solve' the constraints using `smsg` or `smsd`.
For convenience, the Python module has the capacity (on by default) to directly call `smsg` and solve the constraints it produces.

!!! note
    In order for SMS to correctly break symmetries, it is important that Boolean variables that represent graph edges are correctly numbered (see below for the numbering scheme).
    This is best ensured by using the graph encoding-building interface provided by PySMS, and we strongly recommend to use it.

On this page we focus on constraints expressible in propositional logic.
More advanced features are covered in [advanced usage](advanced.md).
This includes generation of graphs with forbidden (induced) subgraphs, constraints specified in custom C++ code, and a technique called [_co-certificate learning_](https://www.ijcai.org/proceedings/2023/216), to solve declarative co-NP complete constraints, such as generating graphs that are _not_ 3-colorable.

## Quick Start: PySMS

If you followed the installation process, the package `pysms` is installed and available for import.

!!! note
	The installed importable module is `pysms`, but the pip package is called `sms-graph-builder`, and it can be removed with `pip uninstall sms-graph-builder` as usual.

The main functionality is located in the `pysms.graph_builder` module, which can be imported for use in scripts (next section), or called directly from the command line, which is the easiest way to start using PySMS.

For example, the following command produces all graphs up to isomorphism with 6 vertices, minimum degree at least 3, and at most 10 edges.

```bash
python -m pysms.graph_builder --vertices 6 --delta-low 3 --num-edges-upp 10 --all-graphs
``` 

The found graphs are printed to standard output as Python lists of edges.

In order to only print the constraints (as a CNF formula in DIMACS) without solving, additionally pass `--no-solve`.
In order to get the results in the `graph6` format (used by [Nauty](https://pallini.di.uniroma1.it/)), pass `--graph-format graph6`.
(Note that this can slow down the computation significantly.)
Below we provide a number of examples to illustrate the functionality.

The most important options to `pysms.graph_builder` are as follows:

- `-v|--vertices n` : search for a graph with `n` vertices;
- `-a|--all-graphs` : enumerate all graphs up to isomorphism (without this, the program terminates after finding the first graph that satisfies the constraints);
- `--no-solve` : don't solve, just output the constraints (possibly to save for a later solve);
- `->|--directed` : generate directed graphs (default is undirected);

and some commonly required bounds:

- `-E|--num-edges-upp` : an upper bound on the number of edges;
- `-e|--num-edges-low` : a lower bound on the number of edges;
- `-D|--Delta-upp` : an upper bound on the maximum degree;
- `-d|--delta-low` : a lower bound on the minimum degree.

Any unrecognized arguments will be forwarded to `smsg`/`smsd` when in solving mode (and ignored with `--no-solve`), so any options that can be set for SMS can also be used through PySMS.

!!! note
	The legacy way of forwarding arguments to SMS via `--args-SMS` is now deprecated.

To get a complete list of all available options for the encoding builder, run 

```bash
python -m pysms.graph_builder --help
```

For all options available for `smsg` or `smsd`, run

```bash
smsg --help
```

!!! tip
    SMS relies on a procedure called the _minimality check_ to filter out non-canonical isomorphic copies of graphs.
    This procedure is often fast, but has worst-case exponential running time.
    To avoid getting stuck in hard corner cases, we strongly recommend to use a time limit (cutoff) for the minimality check, by adding `--cutoff 20000`.
    This limits the number of recursive calls in the minimality check, but potentially results in incomplete symmetry breaking.
    The graphs can be filtered afterwards with tools like [Nauty](https://pallini.di.uniroma1.it/)'s `shortg`.


### Custom Encodings

Next, let us see how to create custom encoding scripts using PySMS.

In the first example, we will create a simple script that uses the `GraphEncodingBuilder` class to create constraints which describe graphs with 7 vertices and maximum degree at most 3.

```python
from pysms.graph_builder import GraphEncodingBuilder
builder = GraphEncodingBuilder(7, directed=False)
builder.maxDegree(3)
builder.solve(allGraphs=True)
```

The `builder` object contains the encoding and all necessary metadata such as the number of vertices and how to map edges to propositional variables. 
The vertices are represented by the integers \(0, 1, ..., n-1\).
A clause can be added using `builder.append(clause)`.

The example above uses only built-in constraints.
The second example below shows how to add custom constraints, in this case to enumerate all *triangle-free* graphs with 7 vertices.

```python
from pysms.graph_builder import GraphEncodingBuilder
from itertools import combinations
builder = GraphEncodingBuilder(7, directed=False)
for i,j,k in combinations(builder.V, 3):
    builder.append([-builder.var_edge(i,j), -builder.var_edge(i,k), -builder.var_edge(j,k)])
builder.solve(allGraphs=True)
```

The function `builder.var_edge(i,j)` returns the Boolean variable associated with the edge \(\{i,j\}\).
So, `[-builder.var_edge(i,j), -builder.var_edge(i,k), -builder.var_edge(j,k)]` represents a clause (a disjunction) that ensures that at least one of the named pairs of vertices has no edge between them, or in other words that the vertex triple `i,j,k` doesn't form a triangle.

!!! note
	The verbose triangle-forbidding code from the previous example can be obtained with the builtin function `GraphEncodingBuilder.ckFree(3)`; we use it for the next example.

The third example goes beyond the builtin functions and encodes _maximal_ triangle-free graphs, i.e. such triangle-free graphs where the addition of any edge creates a triangle (in other words, triangle-free graphs with _diameter_ 2; where any two non-neighbors have a common neighbor).
Notice how in this example we need to create and use _auxiliary_ variables: fresh propositional variables other than those that correspond to edges.
The auxiliary variables are created and returned by `CNF_AND`, and encode the existence of a joint neighbor.

```python
from pysms.graph_builder import GraphEncodingBuilder
from itertools import combinations
builder = GraphEncodingBuilder(7, directed=False)
builder.ckFree(3)
for i,j in combinations(builder.V, 2):
	builder.append(
		[-var_edge(i,j)] +\ # if i and j are non-neighbors, then
		[
            builder.CNF_AND([builder.var_edge(i,k), builder.var_edge(j,k)])
            for k in builder.V if k != i and k != j
        ] # there should be a joint neighbor k
	)
	
builder.solve(allGraphs=True)
```

See [PySMS reference](reference.md) for the full list of builtin constraints and other functions.


## Direct Usage

It is possible to use `smsg` and `smsd` directly, though for them to be useful, one typically needs a set of constraints.
As mentioned earlier, it is strongly recommended to use PySMS to prepare the constraints, although it is perfectly reasonable to store generated constraints and solve them later, or on a different machine.
One particularly useful case is when you want to parallelize the solving of a hard problem, in which case you will likely be calling the SMS solvers directly.

The most important arguments of `smsg` and `smsd` are as follows:

- `-v VERTICES` gives the number of vertices (order) of the graph;
- `-b n m` searches for a bipartite graph with \(n+m\) vertices where the partitions are \(\{0, \dots, n-1\}\) and \(\{n, \dots, n+m-1\}\);
- `--all-graphs` to compute all graphs (instead of stopping at the first solution, which is the default behaviour);
- `--cutoff CUTOFF` limits the number of recursive calls in the minimality check, in order to avoid exponential behaviour. Setting this can render symmetry breaking incomplete.
- `--frequency FREQUENCY` used for balancing the proportion of time spent in the minimality check and in the solver itself. For example if the frequency is 5 then the minimality check is only called every 5-th time (on average; whether to call it is decided by a random coin flip with a `1/FREQUENCY` success rate);
- `--dimacs FILE` the file name providing the SAT encoding in DIMACS format.
- `--print-stats` prints statistics of the run at the end, especially the amount of time spent in various propagators.

SMS has a native interface for the two-step parallelization technique called [cube-and-conquer](https://www.cs.utexas.edu/~marijn/publications/cube.pdf).
In the first step the problem is split into a sequence of subproblems, each represented by a partial assignment (called a _cube_).
The behavior of this phase is controlled by the `--assignment-cutoff*` family of parameters, in particular the two below:

- `--assignment-cutoff ACUTOFF` is used for cubing, if at least `ACUTOFF` edge variables are assigned and propagate is called, then a cube representing the partial assignment is generated.
- `--assignment-cutoff-prerun-time TIME` run for this many seconds before starting cubing (applying assignment cutoff).

In the seconds phase, the cubes are loaded from a file with `--cubes FILE`, and a sub-range to be solved can be picked with `--cube2solve begin end`.

For a complete list of all arguments call `smsg --help` and `smsd --help` respectively.

### Edge Variable Numbering

While it is recommended to use PySMS to build encodings, as long as you adhere to the edge numbering rules (explained below), it is perfectly possible to create valid encodings also with other tools, such as [PySAT](https://pysathq.github.io).

For the undirected version the undirected edge variables are assumed to be given row-wise left-to-right from the upper triangle of the adjacency matrix. For the directed version the directed edge variables are assumed to be given row wise without the diagonal. For example, for an undirected graph on n=4 vertices the corresponding edge variables are given by the following literals: e<sub>0,1</sub> &harr; 1, e<sub>0,2</sub> &harr; 2, e<sub>0,3</sub> &harr; 3, e<sub>1,2</sub> &harr; 4, e<sub>1,3</sub> &harr; 5, e<sub>2,3</sub> &harr; 6. This can also be represented by the following matrix (ignoring diagonal entries):

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

## Examples

Generate all graphs on 7 vertices modulo isomorphism.

```bash
smsg -v 7 --all-graphs
```

Generate all non-bipartite graphs on 7 vertices (with chromatic number at least 3) and show the time spent in the minimality check, and for the coloring, which is performed by a separate propagator:

```bash
smsg -v 7 --all-graphs --min-chromatic-number 3 --print-stats
```

## Solvers

SMS currently requires the SAT solver [CaDiCaL](https://github.com/arminbiere/cadical), but optionally also supports [Clingo](https://potassco.org/clingo).
`build-and-install.sh` will look for Clingo, and if finds it, it will build with Clingo support.
Clingo can be activated with `--use-clingo`.
Note that while SMS started with Clingo, the main development focus now goes to Cadical, and Clingo support is at the moment experimental.
