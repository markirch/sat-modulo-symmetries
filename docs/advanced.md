# Advanced Usage

## Adding custom propagators

Adding custom propagators to SMS is relatively straightforward.
We distinguish between propagators only on fully defined graphs and propagators for partially defined graphs which are only called in a certain frequency.

### Partially defined graphs

For writting a propagator, one can use the class `PartiallyDefinedGraphChecker` and implement the virtual function `virtual void checkProperty(const adjacency_matrix_t &matrix)`.
If the function throws an object of `forbidden_graph_t` given as a list of signed edges, then a clause is automatically added excluding this partially defined graph from the search space.

At the end the an `PartiallyDefinedGraphChecker` object must be added to the solver, i.e., to a  `GraphSolver` object. `GraphSolver` contains a list of `PartiallyDefinedGraphChecker` called `partiallyDefinedGraphCheckers`, where it can simply be added.

#### Example

We give a very simple example, for providing a propagator only keeping triangular free graphs:

```
class TriangleFreeChecker : public PartiallyDefinedGraphChecker
{
    void checkProperty(const adjacency_matrix_t &matrix)
    {
        int n = matrix.size();
        for (int v1 = 0; v1 < n; v1++)
            for (int v2 = v1 + 1; v2 < n; v2++)
                for (int v3 = v2 + 1; v3 < n; v3++)
                {
                    if (matrix[v1][v2] == truth_value_true && matrix[v1][v3] == truth_value_true && matrix[v2][v3] == truth_value_true)
                    {
                        forbidden_graph_t g;
                        g.push_back({truth_value_true, {v1,v2}});
                        g.push_back({truth_value_true, {v1,v3}});
                        g.push_back({truth_value_true, {v2,v3}});
                        throw g;
                    } 
                }
    }
};

```


### Fully defined graphs

Similar, for fully defined graphs, we have the class `FullyDefinedGraphChecker`. The main difference is that the check is only called for fully defined graphs.



## Forbidden Subgraphs

SMS has the optional feature to generate graphs that do not contain any of a set of input forbidden subgraphs.
To enable it, you need to install and build the [Glasgow Subgraph Solver](https://github.com/ciaranm/glasgow-subgraph-solver).
Run the script `pull-and-make-glasgow-solver.sh` to obtain and build the Glasgow Subgraph Solver, then configure with the flag `-DGLASGOW`.

To build with the Glasgow Subgraph Solver, configure with `cmake -Bbuild -H. -DGLASGOW=1` instead (any value for the GLASGOW variable will work).

The forbidden subgraphs are provided to `smsg` by a file. The file is given to the program by the argument `--forbiddenSubgraphs FILE`, where the file contains a forbidden graph in each line. The graphs must be given in an easy to parse format easiest explained on an example

### Example for forbidden subgraphs

```
0 1 1 2 0 2
0 1 1 2 2 3 0 3
```

Two vertices form an edge, i.e., the first line represents the edges {0,1}, {1,2}, {0,2}. So, if we provide this file to the solver, only graphs without a triangle and a four-cycle are generated. Note that the graphs are not checked whether the are induced.

