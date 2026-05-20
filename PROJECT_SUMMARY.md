# SMS Project Summary

## Purpose

SAT Modulo Symmetries (SMS) generates graphs subject to SAT-encoded constraints while avoiding isomorphic duplicates during search. The solver combines a CaDiCaL SAT backend with graph-specific propagators that inspect partial assignments and add blocking clauses.

## Main Layers

- `src/`: C++ solver core, graph encoding, symmetry breaking, and graph property propagators.
- `pysms/`: Python helpers for building DIMACS or QCIR encodings and invoking the solver.
- `encodings/`: Small public example encodings built on top of `pysms`.
- `docs/`: Documentation source for the public project site.
- `cadical_sms/`: CaDiCaL submodule used as the SAT backend.

## Core C++ Flow

1. `src/main.cpp` parses CLI options through `src/options.cpp`.
2. `GraphSolver` in `src/sms.hpp` and `src/sms.cpp` owns the SAT solver, graph handler, and registered graph checkers.
3. `GraphHandler` in `src/graph.hpp` and `src/graph.cpp` maps graph edges to SAT variables and translates assignments to graphs.
4. During search, CaDiCaL callbacks update the current partial assignment.
5. Registered checkers in `src/graphChecker.hpp` inspect partial or complete graphs.
6. Violations are reported by throwing forbidden graph fragments or clauses, which `GraphSolver` adds back to the SAT solver.

## Symmetry Breaking

The symmetry-breaking logic is centered in:

- `src/minimalityCheck.hpp`
- `src/minimalityCheck.cpp`
- `src/minimalityCheckCommon.cpp`
- `src/minimalityCheck_dir.cpp`
- `src/minimalityCheckColex.cpp`

The default minimality check compares graphs using the row-wise ordering. The optional `--colex-ordering` mode uses column-wise/colex ordering. `--sym-break-clauses FILE` can write generated symmetry-breaking clauses in JSON format for inspection or independent checking.

## Solver Options And Features

Important user-facing options are registered in `src/options.cpp`. Common workflows include:

- graph enumeration with `--vertices` and `--all-graphs`
- disabling SMS with `--no-SMS`
- adding DIMACS constraints with `--dimacs`
- enabling property propagators such as connectivity, planarity, forbidden subgraphs, and chromatic number
- creating cubes or simplified formulas for downstream solving
- producing LRAT proofs with `--lrat-output` for UNSAT results

## Python Encoding Layer

`pysms/graph_builder.py` is the main Python entry point for constructing graph SAT encodings. It allocates graph variables, adds clauses for graph properties, writes DIMACS, and can run `smsg`. `pysms/qcir_graph_builder.py` provides similar support for QCIR/QBF workflows.

Useful Python modules include:

- `pysms/solver.py`: Python wrapper around solver execution.
- `pysms/decoder.py`: Converts solver output back to graph edges.
- `pysms/counters.py`: Cardinality/counting helpers used by encodings.

## Build Files

- `build-and-install.sh`: convenience build/install script.
- `CMakeLists.txt` and `src/CMakeLists.txt`: CMake build configuration.
- `pyproject.toml`: Python package metadata.

## Where New Work Should Go

Use this placement rule to keep the public codebase easy to navigate:

- New reusable graph constraints or encoding helpers belong in `pysms/`.
- New example encodings belong in `encodings/`.
- New solver options, propagators, or symmetry-breaking behavior belong in `src/`.
- Detailed user-facing documentation belongs in `docs/`.
- Generated data, solver logs, and one-off experiment outputs should stay outside the repository unless they are intentionally small examples.

Keep `pysms/graph_builder.py` focused on broadly reusable encoding infrastructure. For a specific graph problem, prefer a separate encoding entry point that uses the shared builder API.

## Reproducible Local Runs

When testing branch-specific solver changes, prefer the `smsg` binary and `pysms` module from the current checkout over older installed copies. This avoids accidentally mixing a new encoding with an old solver binary, or the other way around.

## Good Starting Points

For solver behavior, start with `src/sms.cpp`, `src/options.cpp`, and the minimality checker files. For modeling new graph problems, start with `pysms/graph_builder.py` and the examples in `encodings/`.
