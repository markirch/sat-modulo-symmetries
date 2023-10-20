# PySMS Reference

This is a full reference of the Python library PySMS that provides an interface for building encodings of graph problems recognizable by the SMS solvers `smsg` and `smsd`.
The `GraphEncodingBuilder` class provides the interface to build SMS-ready encodings of graph problems, most importantly through `var_edge(u, v)`, which maps graph edges to their propositional variables.
In addition to this, `GraphEncodingBuilder` contains a ready-made library of common constraints, including cardinality constraints.
The module can also be executed from the command line in order to generate or directly solve a formula built from the predifined constraints, which can be passed via command-line parameters.

::: pysms.graph_builder
