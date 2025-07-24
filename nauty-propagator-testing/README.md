# Nauty Propagator Testing

This folder contains testing files and results for the VF2→Nauty replacement in SMS automorphism propagators.

## Background

We replaced the VF2 algorithm with Nauty (dreadnauty) for automorphism counting in SMS propagators:
- **FDG (Fully Defined Graphs)**: Direct Nauty automorphism counting  
- **PDG (Partially Defined Graphs)**: Enumeration of 2^k FDG assignments with Nauty counting

## Test Instance Generation

### CNF File Creation
The test CNF file `cubic_connected_girth4_n12.cnf` was generated using pySMS:

```bash
PYTHONPATH=/Users/szeider/git/discovery/SMS python3 pysms/graph_builder.py \
  --vertices 12 \
  --delta-low 3 --Delta-upp 3 \
  --girth 4 \
  --connectivity-low 1 \
  --cnf-file cubic_connected_girth4_n12.cnf \
  --no-solve
```

**Constraints:**
- **Cubic**: Each vertex has degree exactly 3 (`--delta-low 3 --Delta-upp 3`)
- **Connected**: Graph is connected (`--connectivity-low 1`) 
- **Girth ≥ 4**: No triangles, minimum cycle length 4 (`--girth 4`)
- **Vertices**: n = 12

## Test Results

### Test 1: Baseline (No Automorphism Filtering)
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n12.cnf --vertices 12 --all-graphs --hide-graphs --cutoff 0
```

**Result:** 22 graphs found  
**Expected:** 22 cubic connected girth≥4 graphs for n=12 ✅

### Test 2: FDG Automorphism Filtering  
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n12.cnf --vertices 12 --all-graphs --hide-graphs --cutoff 0 \
  --autcount-frequency 0 --min-automorphisms 64
```

**Result:** 1 graph found  
**Expected:** 1 graph with ≥64 automorphisms (21 graphs filtered out) ✅

**Statistics:**
- AutcountChecker: 22 calls, 21 added clauses (blocking clauses for rejected graphs)
- Only FDG testing (--autcount-frequency 0)
- Perfect filtering behavior

### Test 3: PDG Propagator (Full Integration) - n=12
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n12.cnf --vertices 12 --all-graphs --hide-graphs --cutoff 0 \
  --autcount-frequency 20 --min-automorphisms 64
```

**Result:** 1 graph found  
**Expected:** 1 graph with ≥64 automorphisms ✅

**Statistics:**
- **AutcountPDGChecker**: 83 calls, 47 added clauses (PDG pruning during search)
- **AutcountChecker**: 1 call, 0 added clauses (final FDG check)  
- **MinimalityChecker**: 231 calls (vs 236 in FDG-only mode)
- **Runtime**: 0.18 seconds for PDG checking
- **Search space pruning**: PDG propagator successfully pruned search space

## n=14 Testing Results

### Test 4: Baseline n=14 (No Automorphism Filtering)
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n14.cnf --vertices 14 --all-graphs --hide-graphs --cutoff 0
```

**Result:** 110 graphs found  
**Expected:** 110 cubic connected girth≥4 graphs for n=14 ✅

### Test 5: FDG Automorphism Filtering n=14
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n14.cnf --vertices 14 --all-graphs --hide-graphs --cutoff 0 \
  --autcount-frequency 0 --min-automorphisms 336
```

**Result:** 1 graph found  
**Expected:** 1 graph with ≥336 automorphisms (109 graphs filtered out) ✅

**Statistics:**
- AutcountChecker: 110 calls, 109 added clauses (blocking clauses for rejected graphs)
- Runtime: 0.16 seconds total
- Perfect FDG filtering behavior

### Test 6: PDG Propagator n=14 (Full Integration)
```bash
./build/src/smsg --dimacs cubic_connected_girth4_n14.cnf --vertices 14 --all-graphs --hide-graphs --cutoff 0 \
  --autcount-frequency 20 --min-automorphisms 336
```

**Result:** 1 graph found  
**Expected:** 1 graph with ≥336 automorphisms ✅

**Statistics:**
- **AutcountPDGChecker**: 349 calls, 189 added clauses (PDG pruning during search)
- **AutcountChecker**: 1 call, 0 added clauses (final FDG check)  
- **MinimalityChecker**: 896 calls (vs 915 in FDG-only mode)
- **Runtime**: 0.38 seconds for PDG checking
- **Search space pruning**: Significant pruning with 189 blocked partial graphs

## Implementation Status

### ✅ Completed
- **FDG Nauty Integration**: Direct automorphism counting for fully defined graphs
- **PDG Nauty Integration**: Enumeration of 2^k FDG assignments with Nauty counting
- **FDG Filtering**: Correctly filters graphs based on automorphism count threshold  
- **PDG Search Pruning**: Successfully prunes search space during SMS execution
- **SMS Integration**: Proper integration with both AutcountChecker and AutcountPDGChecker propagators
- **Directed Graph Support**: Fixed Nauty API usage for directed vs undirected graphs
- **Performance Safeguards**: 15-edge limit for PDG enumeration in AutcountPDGChecker

### 🔄 Next Steps  
- **Performance Validation**: Compare PDG enumeration performance vs original VF2
- **Code Cleanup**: Remove dead VF2/igraph code and dependencies
- **Extended Testing**: Test with larger instances and different thresholds

## Key Findings

1. **Complete VF2→Nauty replacement successful** 
   - FDG: Direct Nauty counting works perfectly
   - PDG: 2^k enumeration with Nauty counting works correctly

2. **Search space pruning effective**
   - PDG propagator identified 47 partial graphs that cannot achieve ≥64 automorphisms
   - Reduced MinimalityChecker calls from 236 to 231

3. **Performance acceptable**  
   - FDG testing: ~0.03 seconds
   - PDG testing: ~0.18 seconds (with enumeration overhead)
   - 15-edge limit prevents exponential blowup

4. **SMS integration seamless**
   - Both propagators work correctly with SMS framework
   - Proper clause generation and conflict handling

## Technical Details

**Files Modified:**
- `src/graphPropagators/automorphism_counter.cpp` - Complete Nauty integration
  - `countFDGAutomorphismsNauty()` - Direct FDG counting
  - `nautyPDGThresholdCheck()` - PDG enumeration with early termination
  - `hasAtLeastKAutomorphisms()` - Complete VF2 replacement
  - `countExactAutomorphisms()` - Complete VF2 replacement
- `src/graphPropagators/autcountChecker.cpp` - 15-edge limit in AutcountPDGChecker  
- `src/CMakeLists.txt` - Nauty library integration

**Command Line Options:**
- `--autcount-frequency 0`: FDG-only testing (final graph check)
- `--autcount-frequency N`: PDG testing every N propagator calls  
- `--min-automorphisms K`: Minimum required automorphisms threshold

**Performance Characteristics:**
- **FDG Direct Counting**: O(V!) worst case, very fast in practice
- **PDG Enumeration**: O(2^U) where U = undefined edges, limited to U ≤ 15
- **Early Termination**: Both methods support early termination for threshold checking