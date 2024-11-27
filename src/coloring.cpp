#include <assert.h>
#include "coloring.h"

typedef int color;
typedef int vertex;

typedef struct VCPair
{
	vertex v;
	color c;
} VCPair;

int num_colorings = 0;

// colors the first non-colored vertex with the first available one. Dummy implementation
int simpleColoringAlg(const adjacency_matrix_t &adjacencyMatrix, int maxColors, int nVertices, coloring_t &coloring, int vertex, int highestUsedColor)
{

	if (vertex == nVertices)
		return 1; // all vertices colored

	int highestColorToTry = highestUsedColor < maxColors ? highestUsedColor + 1 : maxColors;
	std::vector<bool> available(highestColorToTry, true);
	for (int i = 0; i < highestColorToTry; i++)
		available[i] = 1;

	for (int u = 0; u < vertex; u++)
		if (adjacencyMatrix[u][vertex])
			available[coloring[u]] = 0;

	for (int i = 0; i < highestColorToTry; i++)
		if (available[i])
		{
			coloring[vertex] = i;
			if (simpleColoringAlg(adjacencyMatrix, maxColors, nVertices, coloring, vertex + 1, i < highestUsedColor ? highestUsedColor : highestColorToTry))
				return 1;
		}
	return 0;
}

// TODO it seems like there was a problem with variable-length arrays, so for now let's do fixed-size
// i.e. we can support up to 32 vertices and 8 colors (we can always increase these constants if necessary)
// this was in C, now we're in C++, so we could also use vectors...
#define VMAX 78
#define CMAX 16

int simpleHyperColoring(const adjacency_matrix_t &adjacencyMatrix, int nVertices[2], coloring_t &coloring, int vertex)
{
	if (vertex == 0)
	{
		coloring[0] = 0;
		return simpleHyperColoring(adjacencyMatrix, nVertices, coloring, 1);
	}
	if (vertex == nVertices[0])
	{
		return 1; // all vertices colored
	}

	int val[2] = {0, 1};
	int can[2] = {1, 1};

	for (int e = nVertices[0]; e < nVertices[0] + nVertices[1] && (can[0] || can[1]); e++)
	{
		if (adjacencyMatrix[vertex][e] == truth_value_true)
		{
			int e_forces[2] = {1, 1};
			for (int v = 0; v < nVertices[0]; v++)
			{
				if (adjacencyMatrix[v][e] == truth_value_true)
				{
					if (v < vertex)
					{
						e_forces[coloring[v]] = 0;
					}
					else if (v > vertex)
					{
						e_forces[0] = e_forces[1] = 0;
						break;
					}
				}
			}
			if (e_forces[0])
				can[1] = 0;
			if (e_forces[1])
				can[0] = 0;
		}
	}

	// TODO some fancy heuristic for the color order
	/*uint32_t cost[2] = {0, 0};

	for (int v = 0; v < vertex; v++) {
	  if (coloring[v] == 0) {
		for (int w = v + 1; w < vertex; w++) {
		  if (coloring[w] == 0) {
			cost[0] += triangle_stats[v][w][vertex];
		  }
		}
	  } else {
		cost[1] += edge_stats[v][vertex];
	  }
	}

	if (cost[1] < cost[0] || (cost[1] == cost[0] && vertex * 3 <= nVertices)) {
	  std::swap(can[0], can[1]);
	  std::swap(val[0], val[1]);
	}*/

	for (int i = 0; i < 2; i++)
	{
		if (can[i])
		{
			coloring[vertex] = val[i];
			if (simpleHyperColoring(adjacencyMatrix, nVertices, coloring, vertex + 1))
				return 1;
		}
	}

	return 0;
}

#include "cadical.hpp"
/**
 * @brief Get a coloring of the graph with a SAT solver if the graph is colorable
 *
 * @param adjacencyMatrix
 * @param coloring
 * @param maxColors
 * @param clique Possible to add a clique at the beginning where the colors should be fixed
 * @return int
 */
bool getColoringSAT(const adjacency_matrix_t &adjacencyMatrix, coloring_t &coloring, int maxColors, vector<int> clique)
{
	int var_counter = 1;
	vector<vector<int>> colors;
	int nVertices = adjacencyMatrix.size();

	auto solver = CaDiCaL::Solver();

	for (int i = 0; i < nVertices; i++)
	{
		colors.push_back(vector<int>());
		for (int l = 0; l < maxColors; l++)
			colors[i].push_back(var_counter++);

		// at least one color
		for (auto l : colors[i])
			solver.add(l);
		solver.add(0);
	}

	for (int i = 0; i < nVertices; i++)
		for (int j = i + 1; j < nVertices; j++)
		{
			if (adjacencyMatrix[i][j] == truth_value_true)
				for (int c = 0; c < maxColors; c++)
				{
					solver.add(-colors[i][c]);
					solver.add(-colors[j][c]);
					solver.add(0);
				}
		}

	if ((int) clique.size() > maxColors)
		return false; // uncolorable if there is a large clique
	for (int i = 0; i < (int) clique.size(); i++)
	{
		solver.add(colors[clique[i]][i]);
		solver.add(0);
	}

	// some symmetry breaking smallest available color; i.e.; for each smaller color than the selected one: adjacent to smaller vertex with this color or adjacent to a vertex in the clique with this color
	// TODO

	auto res = solver.solve();
	if (res == 10)
	{
		// extract coloring
		for (int i = 0; i < nVertices; i++)
			for (int c = 0; c < maxColors; c++)
			{
				if (solver.val(colors[i][c]) > 0)
				{
					coloring[i] = c;
					break;
				}
			}
	}
	return res == 10;
}

int simpleColoring010(const adjacency_matrix_t &adjacencyMatrix, int nVertices, coloring_t &coloring, int vertex, vector<vector<vector<int>>> &triangle_stats, vector<vector<int>> &edge_stats)
{
	if (vertex == nVertices)
	{
		return 1; // all vertices colored
	}

	int val[2] = {0, 1};
	int can[2] = {1, 1};

	for (int u = 0; u < vertex && (can[0] + can[1] > 0); u++)
	{
		if (adjacencyMatrix[u][vertex])
		{
			if (coloring[u] == 1)
			{
				can[1] = 0;
			}
			else
			{
				for (int v = u + 1; v < vertex; v++)
				{
					if (adjacencyMatrix[u][v] && adjacencyMatrix[v][vertex] && coloring[v] == 0)
					{
						can[0] = 0;
					}
				}
			}
		}
	}

	/* count the number of vertices colored 1
	int t = 0;
	for (int v = 0; v < vertex; v++) {
		t += coloring[v];
	}*/

	// an experimentally found heuristic that tends to produce a small number of colorings overall
	/*if (vertex <= nVertices / 3 + 2 * (nVertices % 5 - nVertices % 3)) {
		std::swap(can[0], can[1]);
		std::swap(val[0], val[1]);
	}*/

	uint32_t cost[2] = {0, 0};

	for (int v = 0; v < vertex; v++)
	{
		if (coloring[v] == 0)
		{
			for (int w = v + 1; w < vertex; w++)
			{
				if (coloring[w] == 0)
				{
					cost[0] += triangle_stats[v][w][vertex];
				}
			}
		}
		else
		{
			cost[1] += edge_stats[v][vertex];
		}
	}

	if (cost[1] < cost[0] || (cost[1] == cost[0] && vertex * 3 <= nVertices))
	{
		std::swap(can[0], can[1]);
		std::swap(val[0], val[1]);
	}

	for (int i = 0; i < 2; i++)
	{
		if (can[i])
		{
			coloring[vertex] = val[i];
			if (simpleColoring010(adjacencyMatrix, nVertices, coloring, vertex + 1, triangle_stats, edge_stats))
				return 1;
		}
	}

	return 0;
}

int getHyperColoring(int nVertices[2], const adjacency_matrix_t &adjacencyMatrix, coloring_t &coloring)
{
	int colorable = simpleHyperColoring(adjacencyMatrix, nVertices, coloring, 0);
	num_colorings += colorable;
	return colorable;
}

vector<lit_t> getHyperColoringClause(const coloring_t &coloring, int nVertices[2], const adjacency_matrix_t &matrix, const vector<vector<lit_t>> &E)
{
	vector<lit_t> clause;
	for (int e = nVertices[0]; e < nVertices[0] + nVertices[1]; e++)
	{
		int witness[2] = {-1, -1};
		for (int v = 0; v < nVertices[0]; v++)
		{
			if (matrix[v][e] == truth_value_true)
			{
				witness[coloring[v]] = v;
				if (witness[0] != -1 && witness[1] != -1)
				{
					break;
				}
			}
		}
		if (witness[0] == -1)
		{
			printf("error 0\n");
		}
		clause.push_back(-E[witness[0]][e]);
		if (witness[1] == -1)
		{
			printf("error 1\n");
		}
		clause.push_back(-E[witness[1]][e]);
	}
	return clause;
}

vector<vector<int>> extractColorClasses(const coloring_t &coloring, int nVertices)
{
	vector<vector<int>> color_class(2);
	for (int v = 0; v < nVertices; v++)
	{
		if (coloring[v] >= (int) color_class.size())
		{
			color_class.resize(coloring[v] + 1);
		}
		color_class[coloring[v]].push_back(v);
	}
	return color_class;
}

// returns a Tseitin-converted circuit that say the coloring should not work
// the circuit says that at least one edge should be monochromatic under the coloring
vector<vector<lit_t>> getHyperColoringCircuit(const coloring_t &coloring, int nVertices[2], const vector<vector<lit_t>> &E, int &next_var)
{
	vector<vector<int>> color_class = extractColorClasses(coloring, nVertices[0]);

	vector<vector<lit_t>> circuit(1);
	for (int e = nVertices[0]; e < nVertices[0] + nVertices[1]; e++)
	{
		for (int color = 0; color < 2; color++)
		{
			int monochromatic = next_var++;
			circuit.push_back({+monochromatic});
			for (int v : color_class[1 - color])
			{
				circuit.back().push_back(+E[v][e]);
			}
			for (int v : color_class[1 - color])
			{
				circuit.push_back({-monochromatic, -E[v][e]});
			}
			circuit.front().push_back(monochromatic);
		}
	}
	return circuit;
}

int get010Coloring(int nVertices, const adjacency_matrix_t &adjacencyMatrix, coloring_t &coloring, vector<vector<vector<int>>> triangle_stats, vector<vector<int>> edge_stats)
{
	int colorable = simpleColoring010(adjacencyMatrix, nVertices, coloring, 0, triangle_stats, edge_stats);
	num_colorings += colorable;
	return colorable;
}

vector<lit_t> get010ColoringClause(const coloring_t &coloring, int nVertices, const vector<vector<lit_t>> &E, const vector<vector<vector<lit_t>>> &T)
{
	vector<lit_t> clause;
	for (int u = 0; u < nVertices; u++)
	{
		if (coloring[u] == 1)
		{
			for (int v = u + 1; v < nVertices; v++)
				if (coloring[v] == 1)
					clause.push_back(E[u][v]);
		}
		else
		{
			for (int v = u + 1; v < nVertices; v++)
				if (coloring[v] == 0)
					for (int w = v + 1; w < nVertices; w++)
						if (coloring[w] == 0)
							clause.push_back(T[u][v][w]);
		}
	}
	return clause;
}

int coloringDPLL(int vertex_heuristic, int color_heuristic, const adjacency_matrix_t &adjacencyMatrix, int maxColors, int nVertices, coloring_t &coloring)
{
	unsigned int colorUsageCount[CMAX];
	vertex trace[VMAX];
	int isDecision[VMAX];
	int numColored = 0;
	unsigned int decLevel = 0;
	unsigned int vertexColorBlockingLevel[VMAX][CMAX];
	unsigned int nextColor[VMAX];
	unsigned int numAvailable[VMAX];
	color availableColors[VMAX][CMAX];
	VCPair propagationStack[VMAX];
	unsigned int propStackSize = 0;

	// TODO implement backjumping (probably need to keep a reason vertex for every blocked color)

	for (vertex v = 0; v < nVertices; v++)
	{
		coloring[v] = -1;
		trace[v] = -1;
		nextColor[v] = 0;
		for (color c = 0; c < maxColors; c++)
		{
			vertexColorBlockingLevel[v][c] = -1;
		}
	}
	for (color c = 0; c < maxColors; c++)
	{
		colorUsageCount[c] = 0;
	}

	propagationStack[propStackSize++] = {0, 0};
	colorUsageCount[0]++;

	// symmetry breaking: color vertex 0 with 0, and 0's first neighbor with 1 (will fail with maxColors=1)
	// TODO experiment more with vertices and colors to be picked here
	/* doesn't seem to work very well, so turned off for now
	vertex v0 = 0, v1 = 0;
	// v0 := vertex with max degree
	// v1 := v0's neighbor with max degree
	unsigned int D = 0;
	for (vertex v = 0; v < nVertices; v++) {
		unsigned int deg = 0;
		for (vertex w = 0; w < nVertices; w++) {
			deg += adjacencyMatrix[v][w];
		}
		if (deg > D) {
			v0 = v;
			D = deg;
		}
	}
	D = 0;
	for (vertex v = 0; v < nVertices; v++) {
		if (adjacencyMatrix[v][v0]) {
			unsigned int deg = 0;
			for (vertex w = 0; w < nVertices; w++) {
				deg += adjacencyMatrix[v][w];
			}
			if (deg > D) {
				v1 = v;
				D = deg;
			}
		}
	}
	color c0 = 0, c1 = 1;

	//propagationStack[propStackSize++] = {v0, c0};
	//propagationStack[propStackSize++] = {v1, c1};
	*/

	while (numColored < nVertices)
	{
		// unitPropagate();
		while (propStackSize > 0)
		{
			VCPair vc = propagationStack[--propStackSize]; // pop
			// colorVertex(vc.v, vc.c, coloring, trace, &numColored);
			// fprintf(stderr, "coloring %d with color %d\n", vc.v, vc.c);
			coloring[vc.v] = vc.c;
			trace[numColored++] = vc.v;
			// colorUsageCount[vc.c]++;
			for (vertex w = 0; w < nVertices; w++)
			{
				// fprintf(stderr, "look, I just colored %d with %d, ok? Now I am looking at %d (neighbour=%d, coloring=%d, VCBL=%u)\n", vc.v, vc.c, w, adjacencyMatrix[w][vc.v], coloring[w], vertexColorBlockingLevel[w][vc.v]);
				if (
					adjacencyMatrix[w][vc.v] &&
					coloring[w] == -1 &&
					vertexColorBlockingLevel[w][vc.c] > decLevel)
				{
					// fprintf(stderr, "disabling %d for the neighbour %d\n", vc.c, w);
					//  color vc.c is unavailable for vertex w starting with decision level decLevel
					vertexColorBlockingLevel[w][vc.c] = decLevel;

					// find the lowest and highest available color for w (to determine if w is unit or conflicting
					color l = 0, h = maxColors - 1;
					while (l < maxColors && vertexColorBlockingLevel[w][l] <= decLevel)
					{
						++l;
					}
					if (l == maxColors)
					{
						// conflict, backtrack
						// fprintf(stderr, "conflict at decision level %u\n", decLevel);
						if (decLevel == 0)
						{
							// UNSAT
							return 0;
						}
						// delete propagation stack
						for (unsigned int i = 0; i < propStackSize; i++)
						{
							colorUsageCount[propagationStack[i].c]--;
						}
						propStackSize = 0;

						// pop all propagated things from the trace
						while (!isDecision[trace[--numColored]])
						{
							vertex p = trace[numColored];
							color c = coloring[p];
							coloring[p] = -1;
							colorUsageCount[c]--;

							// unblock c for p's neighbours for which it was blocked at this decision level
							for (vertex w = 0; w < nVertices; w++)
							{
								if (
									adjacencyMatrix[w][p] == 1 &&
									coloring[w] == -1 &&
									vertexColorBlockingLevel[w][c] == decLevel)
								{
									vertexColorBlockingLevel[w][c] = -1;
								}
							}
						}

						// we are now at the last decision
						vertex lastDec = trace[numColored];
						color lastC = coloring[lastDec];
						color nextC = availableColors[lastDec][nextColor[lastDec]]; // this is guaranteed to be a valid color
						coloring[lastDec] = -1;
						colorUsageCount[lastC]--;
						assert(nextC < maxColors);
						assert(nextC != lastC);
						nextColor[lastDec]++;
						// fprintf(stderr, "last decision was vertex %d with color %d\n", lastDec, lastC);
						// fprintf(stderr, "will try %d with color %d next\n", lastDec, nextC);
						//  if the next color is the last available, change this decision to propagation
						unsigned int newDecLevel = decLevel;
						if (nextColor[lastDec] == numAvailable[lastDec])
						{
							newDecLevel--;
							// fprintf(stderr, "this will be propagated\n");
						}
						else
						{
							// fprintf(stderr, "this will be another decision\n");
						}

						// unblock lastC for lastDec's neighbours for which it was blocked at this decision level
						for (vertex w = 0; w < nVertices; w++)
						{
							if (
								adjacencyMatrix[w][lastDec] == 1 &&
								coloring[w] == -1 &&
								vertexColorBlockingLevel[w][lastC] == decLevel)
							{
								vertexColorBlockingLevel[w][lastC] = -1;
							}
						}
						if (decLevel > newDecLevel)
						{
							isDecision[lastDec] = 0;
							decLevel--;
						}
						propagationStack[propStackSize++] = {lastDec, nextC};
						colorUsageCount[nextC]++;
						break; // from the for cycle that updates colors for neighbours of an assigned vertex (and into the propagation loop)
					}
					else
					{
						// find the highest color that is available and, if l is fresh, not fresh
						if (colorUsageCount[l] == 0)
						{
							while (l < h && (vertexColorBlockingLevel[w][h] <= decLevel || colorUsageCount[h] == 0))
							{
								--h;
							}
						}
						else
						{
							while (l < h && (vertexColorBlockingLevel[w][h] <= decLevel))
							{
								--h;
							}
						}
						if (l == h)
						{
							// unit
							// fprintf(stderr, "the vertex %d is now unit with color %d\n", w, l);
							propagationStack[propStackSize++] = {w, l};
							colorUsageCount[l]++;
							isDecision[w] = 0;
						}
					}
				}
			}
		}

		/* propagation reached fixpoint
		 * find a vertex with the smallest domain
		 * and make a decision */
		if (numColored == nVertices)
		{
			break;
		}
		// fprintf(stderr, "--- very good, we have already colored %u/%u vertices\n", numColored, nVertices);
		/* calculate the vertices with
		 *   minumum number of available colors (mdv = min domain vertex)
		 *   maximum number of available colors (Mdv = max domain vertex)
		 *   minimum number of uncolored neighbors (munv = min uncolored neighbors vertex)
		 *   maximum number of uncolored neighbors (Munv = min uncolored neighbors vertex)
		 */
		vertex mdv = -1, Mdv = -1, munv = -1, Munv = -1, rv = -1;
		unsigned int minDomainSize = maxColors + 1, maxDomainSize = 1;
		unsigned int minUncoloredNeighbors = nVertices, maxUncoloredNeighbors = 0;
		unsigned int numUncoloredSeen = 0;
		for (vertex v = 0; v < nVertices; v++)
		{
			if (coloring[v] == -1)
			{
				unsigned int domainSize = 0;
				for (color c = 0; c < maxColors; c++)
				{
					// fprintf(stderr, "random log v=%u c=%d\n", v, c);
					if (vertexColorBlockingLevel[v][c] > decLevel)
					{
						domainSize++;
					}
				}
				unsigned int uncoloredNeighbors = 0;
				for (vertex w = 0; w < nVertices; w++)
				{
					if (adjacencyMatrix[v][w] && coloring[w] == -1)
					{
						uncoloredNeighbors++;
					}
				}
				if (domainSize < minDomainSize)
				{
					mdv = v;
					minDomainSize = domainSize;
				}
				if (domainSize > maxDomainSize)
				{
					Mdv = v;
					maxDomainSize = domainSize;
				}
				if (uncoloredNeighbors < minUncoloredNeighbors)
				{
					munv = v;
					minUncoloredNeighbors = uncoloredNeighbors;
				}
				if (uncoloredNeighbors >= maxUncoloredNeighbors)
				{
					Munv = v;
					maxUncoloredNeighbors = uncoloredNeighbors;
				}
				// the k-th seen vertex is picked with probability 1/k -> will result in uniform distribution
				if (random() % ++numUncoloredSeen == 0)
				{
					rv = v;
				}
			}
		}
		assert(mdv > -1);
		assert(Mdv > -1);
		assert(munv > -1);
		assert(Munv > -1);
		assert(rv > -1);

		// rationale behind choosing Munv: vertices with few uncolored neighbors won't have an impact, so we can color them later
		vertex nextDec = mdv;

		switch (vertex_heuristic)
		{
		case 0:
			nextDec = mdv;
			break;
		case 1:
			nextDec = Mdv;
			break;
		case 2:
			nextDec = munv;
			break;
		case 3:
			nextDec = Munv;
			break;
		case 4:
			nextDec = rv;
			break;
		}

		// fprintf(stderr, "nextDec=%d\n", nextDec);
		int have_new_color = 0;
		numAvailable[nextDec] = 0;
		for (color c = 0; c < maxColors; c++)
		{
			if (vertexColorBlockingLevel[nextDec][c] > decLevel)
			{
				// symmetry breaking: only consider one unused color
				if (colorUsageCount[c] == 0)
				{
					if (have_new_color)
					{
						continue;
					}
					have_new_color = 1;
				}
				availableColors[nextDec][numAvailable[nextDec]++] = c;
			}
		}

		switch (color_heuristic)
		{
		case 0: // sort in increasing order of usage
			for (unsigned int i = 1; i < numAvailable[nextDec]; i++)
			{
				unsigned int j = i;
				color tc = availableColors[nextDec][i];
				unsigned int tcUse = colorUsageCount[tc];
				while (j > 0 && colorUsageCount[availableColors[nextDec][j - 1]] > tcUse)
				{
					availableColors[nextDec][j] = availableColors[nextDec][j - 1];
					--j;
				}
				availableColors[nextDec][j] = tc;
			}
			break;
		case 1: // sort in decreasing order of usage
			for (unsigned int i = 1; i < numAvailable[nextDec]; i++)
			{
				unsigned int j = i;
				color tc = availableColors[nextDec][i];
				unsigned int tcUse = colorUsageCount[tc];
				while (j > 0 && colorUsageCount[availableColors[nextDec][j - 1]] < tcUse)
				{
					availableColors[nextDec][j] = availableColors[nextDec][j - 1];
					--j;
				}
				availableColors[nextDec][j] = tc;
			}
			break;
		case 2: // shuffle randomly
			for (unsigned int i = numAvailable[nextDec]; i > 0; i--)
			{
				int ridx = random() % i;
				color tc = availableColors[nextDec][i - 1];
				availableColors[nextDec][i - 1] = availableColors[nextDec][ridx];
				availableColors[nextDec][ridx] = tc;
			}
			break;
		}

		// fprintf(stderr, "random log nextDec=%d\n", nextDec);
		nextColor[nextDec] = 1;
		propagationStack[propStackSize++] = {nextDec, availableColors[nextDec][0]};
		colorUsageCount[availableColors[nextDec][0]]++;
		isDecision[nextDec] = 1;
		decLevel++;
	}

	return 1;
}

// TODO rewrite subfunctions
// returns empty coloring if no coloring was found
bool Coloring::getColoring(int nVertices, const adjacency_matrix_t &adjacencyMatrix, coloring_t &coloring, int coloringAlgo, vector<int> clique)
{
	// printf("Start\n");
	// clock_t start = clock();
#ifndef DIRECTED
	if (coloringAlgo == 1)
	{
		return coloringDPLL(3, 1, adjacencyMatrix, maxColors, nVertices, coloring);
	}
	else if (coloringAlgo == 2)
	{
		return getColoringSAT(adjacencyMatrix, coloring, maxColors, clique);
	}
	else
	{
		return simpleColoringAlg(adjacencyMatrix, maxColors, nVertices, coloring, 0, 0);
	}
#else
	// we use colorable if the undirected version is uncolorable
	auto m_cpy = adjacencyMatrix;
	for (int i = 0; i < nVertices; i++)
	{
		for (int j = i + 1; j < nVertices; j++)
		{
			if (adjacencyMatrix[i][j] == truth_value_true || adjacencyMatrix[j][i] == truth_value_true)
			{
				m_cpy[i][j] = m_cpy[j][i] = truth_value_true;
			}
		}
	}

	if (coloringAlgo == 1)
	{
		return coloringDPLL(3, 1, m_cpy, maxColors, nVertices, coloring);
	}
	else if (coloringAlgo == 2)
	{
		return getColoringSAT(m_cpy, coloring, maxColors, clique);
	}
	else
	{
		return simpleColoringAlg(m_cpy, maxColors, nVertices, coloring, 0, 0);
	}
#endif
}

std::vector<edge_t> Coloring::coloring2monochromaticEdges(const coloring_t &coloring)
{
	std::vector<edge_t> edgeList;

	for (size_t i = 0; i < coloring.size(); i++)
		for (size_t j = i + 1; j < coloring.size(); j++)
			if (coloring[i] == coloring[j])
			{
				edgeList.push_back(std::make_pair(i, j));
#ifdef DIRECTED
				edgeList.push_back(std::make_pair(j, i)); // both directions
#endif
			}
	return edgeList;
}
