# Advanced Usage

## Adding custom propagators

Adding custom propagators to SMS is relatively straightforward.
We distinguish between propagators only on fully defined graphs and propagators for partially defined graphs which are only called in a certain frequency.

### Partially defined graphs

For writing a propagator, one can use the class `PartiallyDefinedGraphChecker` and implement the virtual function `virtual void checkProperty(const adjacency_matrix_t &matrix)`.
If the function throws an object of `forbidden_graph_t` given as a list of signed edges, then a clause is automatically added excluding this partially defined graph from the search space.

At the end, a `PartiallyDefinedGraphChecker` object must be added to the solver, i.e., to a  `GraphSolver` object. `GraphSolver` contains a list of `PartiallyDefinedGraphChecker` called `partiallyDefinedGraphCheckers`, where it can simply be added.

#### Example

We give a very simple example, of a propagator that only keeps triangle-free graphs:

```python
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

SMS has the optional feature to generate graphs that do not contain any of a set of input forbidden graphs, either as ordinary or induced subgraphs.
To enable the feature, install and build with the [Glasgow Subgraph Solver](https://github.com/ciaranm/glasgow-subgraph-solver), by running `./build-and-install.sh -s`.

The forbidden subgraphs are provided to `smsg` by a file. The file is given to the program by the argument `--forbidden-subgraphs FILE` (or `--forbidden-induced-subgraphs FILE`, the two can be used together), where the file contains a forbidden graph in each line.
Each line describes a graph as follows: the line must contain an odd number of space-separated non-negative integers; the first one gives the number of vertices, and each following consecutive pair describes one edge.

### Example for forbidden subgraphs

```plaintext
3 0 1 1 2 0 2
3 0 1 1 2 2 3 0 3
```

Each consecutive pair of vertices form an edge, i.e., the first line represents a graph with 3 vertices and the edges {0,1}, {1,2}, {0,2}. So, if we provide this file to the solver with `--forbidden-subgraphs`, only graphs without triangles and four-cycles will be generated.

