# TDA Flood Simulation — Coding Plan & Claude Code Prompts

**Project**: Topological flood simulation, Biesbosch–Noordwaard (AHN DTM)  
**Goal**: Implement Morse-theory-based flood analysis, comparing simple threshold vs. connectivity-based topological simulation

---

## 1. Current notebook diagnosis

The notebook (`vis_v2_data_changed.ipynb`) already has solid data infrastructure. Here's what's done and what's missing:

### Done
- Loads and mosaics 4 AHN `.TIF` tiles using `rasterio`
- Histogram, terrain map, contour visualization
- `find_local_minima_simple` — naive 4-neighbor Python loop (correct but very slow)
- `find_meaningful_minima` — Gaussian smoothing + `h_minima` (better, but orphaned from the simulation)
- `flood_fill` — BFS from a single hardcoded origin pixel

### Critical issues to fix before extending
| Issue | Problem | Fix |
|---|---|---|
| `flood_fill` origin is arbitrary `(2000, 500)` | Not hydrologically meaningful | Derive seed pixels from the data (NaN/low-elevation = river) |
| NaN filled with `nanmin - 1.0` | Creates artificial ultra-low cells that dominate any topological sweep | Fill with a large "barrier" value instead, or keep as masked |
| `threshold` in `flood_fill` means height tolerance above origin, not a rising water level | Conceptually wrong for a flood simulation | Replace with proper sublevel set / priority-queue approach |
| `find_local_minima_simple` uses Python loops | Prohibitively slow on 3000×3000 | Replace with numpy/scipy vectorized approach |
| Persistence filtering (`find_meaningful_minima`) not connected to the simulation | Unused result | Wire persistence threshold into the merge tree and flood sim |

---

## 2. Recommended packages

| Package | Why | Install |
|---|---|---|
| `rasterio`, `numpy`, `scipy`, `skimage`, `matplotlib` | Already there | — |
| **`gudhi`** | CubicalPersistence: directly computes persistence pairs on a 2D grid (birth = basin min elevation, death = saddle elevation). The principled replacement for `h_minima`. Also gives saddle *locations*. | `pip install gudhi` |
| **`scipy.cluster.hierarchy.DisjointSet`** | Built-in union-find. Used to build the merge tree via the sublevel set sweep. No extra install. | — |
| **`heapq`** (stdlib) | Min-heap for priority-queue flood simulation (the hydrologically correct way to grow flooded regions from seeds). | — |
| **`networkx`** | Build and visualize the merge tree as an actual graph/dendrogram | `pip install networkx` |
| `matplotlib.animation.FuncAnimation` | Flood animation (already imported) | — |

### Why not TTK?
TTK is the gold standard for topology in scientific visualization, with full merge tree / contour tree / Morse-Smale complex support. But it requires ParaView integration or a complex C++ build — not Jupyter-friendly. GUDHI + the custom union-find sweep below gives you everything you need for this project and is fully explainable in the report.

### Why GUDHI over just `h_minima`?
`h_minima` finds depressions of depth ≥ h, but it doesn't give you:
- The saddle elevation where two basins merge
- The merge tree structure
- A persistence diagram for principled threshold selection

GUDHI's `CubicalComplex` computes all of this in one call. The persistence pairs `(birth, death)` are exactly `(min_elevation, saddle_elevation)`, and persistence = `death - birth` = basin depth (= the same `h` as in `h_minima`, but now computed for all basins at once).

---

## 3. Architecture overview

```
AHN TIF tiles
     │
     ▼
[Module 1] Data prep
  - Mosaic + crop
  - NaN audit: identify river seed pixels (NaN or elevation < 0m NAP)
  - Fill remaining NaN with max_elevation + 1 (barrier)
     │
     ├──────────────────────────────────────┐
     ▼                                      ▼
[Module 2] GUDHI persistence          [Module 3] Union-Find merge tree
  - CubicalComplex on crop               - Sort all pixels by elevation
  - Persistence pairs (b, d)             - Sweep from low to high
  - Persistence diagram                  - DisjointSet: track components
  - Filter noise by persistence          - Record merge events (saddles)
  - Identify saddle pixel locations      - Output: merge tree graph
     │                                      │
     └───────────────┬──────────────────────┘
                     ▼
           [Module 4] Topological flood sim
             - Seeds: river pixels identified in Module 1
             - Priority queue: grow flooded region by ascending elevation
             - Label each pixel with its source basin
             - Record component merge events
                     │
          ┌──────────┴──────────┐
          ▼                     ▼
  [Module 5] Comparison   [Module 6] Merge tree viz
    Simple threshold vs      Dendrogram plot
    topological flood        Saddle points on terrain
          │
          ▼
  [Module 7] Animation
    FuncAnimation of
    rising water level
```

---

## 4. Module breakdown & Claude Code prompts

Each prompt is self-contained and references the existing notebook. Paste them directly into Claude Code.

---

### MODULE 1 — NaN handling & seed identification

**What it does**: Audits the NaN situation, derives hydrologically meaningful seed pixels (river/water entry points), and prepares a "clean" elevation array for topological analysis.

```
PROMPT 1 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
I have a Jupyter notebook with an AHN elevation array called `crop` 
(shape ~3000×3000, dtype float64, NAP-relative elevation in metres).
Many cells are NaN — these correspond to water bodies (rivers, creeks) 
in the Biesbosch area. I also have `crop_filled` which fills NaN with 
nanmin-1.0, but this is problematic for topological analysis.

Please add a new notebook section called "Data preparation for TDA" 
that does the following:

1. SEED IDENTIFICATION: Find all NaN pixels and all pixels with 
   elevation < 0.0 m NAP. These are the "river seed pixels" — entry 
   points from which flooding originates. 
   - Create a boolean mask called `seed_mask` (True where seeds are).
   - Print how many seed pixels were found and what fraction of the crop they are.
   - Plot the seed pixels overlaid on the terrain (blue on grey terrain).

2. BARRIER FILL: Create `crop_tda` by copying `crop` and replacing NaN 
   with `np.nanmax(crop) + 1.0`. This makes NaN cells into "walls" for 
   the topological analysis (they will never be in a sublevel set).
   - Keep `crop_filled` as-is for backwards compatibility with existing cells.
   - Print min/max/mean of `crop_tda`, confirming no NaN remain.

3. CONNECTIVITY CHECK: Using `scipy.ndimage.label`, check if the seed 
   pixels form one connected component or several. Print the number of 
   seed components and show their bounding boxes.

Variables to produce: `seed_mask` (bool array), `crop_tda` (float array, 
no NaN). Both same shape as `crop`.
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 2 — GUDHI persistence diagram

**What it does**: Computes persistent homology on the elevation grid using GUDHI's `CubicalComplex`. Gives a persistence diagram, shows which basins are noise vs. real, and locates saddle points. This replaces and extends `find_meaningful_minima`.

```
PROMPT 2 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
I have a 2D numpy array `crop_tda` (elevations in metres, no NaN, 
NAP-relative). I want to use GUDHI's CubicalComplex to compute persistent 
homology of the sublevel sets — this will give me the "birth" and "death" 
elevations of every basin in the terrain.

Add a notebook section called "Persistence analysis (GUDHI)" that does:

1. GUDHI SETUP: Install gudhi if not present. Import 
   `gudhi.CubicalComplex`. Create a CubicalComplex from `crop_tda` 
   (pass it as `top_dimensional_cells=crop_tda.flatten()` with 
   `dimensions=list(crop_tda.shape)`).
   Note: GUDHI uses the array as a cubical complex where each cell's 
   value is the elevation. Sublevel set persistence gives H0 pairs 
   (connected components) which are exactly the basin birth/death pairs.

2. COMPUTE PERSISTENCE: Call `.persistence()` on the complex.
   Filter to dimension 0 (connected components). Extract birth/death 
   arrays. Persistence = death - birth = basin depth.
   Print: number of basins found, min/max/mean persistence.

3. PERSISTENCE DIAGRAM: Plot the persistence diagram (birth on x-axis, 
   death on y-axis, diagonal y=x). Color points by persistence value 
   (colormap = 'plasma'). Mark the diagonal clearly.
   Add a slider or text input to show how many basins survive above a 
   given persistence threshold.

4. PERSISTENCE HISTOGRAM: Plot a histogram of persistence values to 
   help choose a noise-filtering threshold. Mark suggested thresholds 
   at 0.1m, 0.25m, 0.5m, 1.0m.

5. BASIN LOCATIONS: For a chosen threshold (default 0.5m), identify 
   the pixel locations of surviving minima. GUDHI gives you birth 
   simplex indices — convert these back to (row, col) coordinates using 
   `np.unravel_index`. Plot these on the terrain.

6. SADDLE LOCATIONS: Similarly identify the death simplex locations 
   (the saddle points at which basins merge). Plot saddles on the terrain 
   with a different marker.

Use `persistence_threshold = 0.5` as the default; make it easy to change.
Output variables: `persistence_pairs` (Nx2 array of birth/death values), 
`basin_locations` (list of (row,col) for minima above threshold), 
`saddle_locations` (list of (row,col) for corresponding saddle points).
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 3 — Union-Find merge tree (sublevel set sweep)

**What it does**: Implements the core Morse theory algorithm from scratch using `scipy.cluster.hierarchy.DisjointSet`. Sweeps elevation levels from low to high, merges components at saddle points, and builds the merge tree. This is the pedagogically important part — it shows *how* the topology changes.

```
PROMPT 3 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
I want to implement a sublevel set sweep to build the merge tree of 
the terrain `crop_tda` from scratch using union-find. This is the 
Morse theory algorithm: as we "raise the water level", new components 
are born at local minima and merge at saddle points.

Add a notebook section called "Merge tree via sublevel set sweep" with:

1. ALGORITHM OVERVIEW (comment block):
   - Sort all pixels by elevation (ascending)
   - Process each pixel in order:
     * The pixel "appears" (potentially births a new component)
     * Check its already-processed 8-neighbors
     * If they belong to different components: this pixel is a saddle — 
       merge those components via union-find
     * If no processed neighbors: this pixel is a local minimum — birth 
       a new component
     * Record birth events as (elevation, pixel_row, pixel_col)
     * Record merge events as (elevation, component_A_birth, component_B_birth)

2. IMPLEMENTATION:
   Use `scipy.cluster.hierarchy.DisjointSet` (available in scipy >= 1.6).
   
   ```python
   from scipy.cluster.hierarchy import DisjointSet
   
   def compute_merge_tree(Z, connectivity=8):
       """
       Z: 2D elevation array (no NaN)
       Returns: births list, merges list, component_label_array
       """
       rows, cols = Z.shape
       flat = Z.flatten()
       order = np.argsort(flat)  # pixels sorted by elevation, low to high
       
       ds = DisjointSet(range(rows * cols))
       processed = np.zeros(rows * cols, dtype=bool)
       
       # birth_elev[pixel_idx] = elevation at which this pixel's component was born
       birth_elev = np.full(rows * cols, np.nan)
       
       births = []   # (elevation, row, col)
       merges = []   # (saddle_elevation, component_root_A, component_root_B)
       
       for idx in order:
           r, c = divmod(idx, cols)
           processed[idx] = True
           
           # find processed neighbors
           neighbor_roots = set()
           for dr, dc in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]:
               nr, nc = r+dr, c+dc
               if 0 <= nr < rows and 0 <= nc < cols:
                   nidx = nr * cols + nc
                   if processed[nidx]:
                       neighbor_roots.add(ds[nidx])
           
           if len(neighbor_roots) == 0:
               # local minimum: new component born
               birth_elev[idx] = flat[idx]
               births.append((flat[idx], r, c))
           else:
               # possibly a saddle: merge components
               roots_list = list(neighbor_roots)
               for other_root in roots_list[1:]:
                   if ds[roots_list[0]] != ds[other_root]:
                       merges.append((
                           flat[idx],
                           birth_elev[ds[roots_list[0]]],
                           birth_elev[ds[other_root]]
                       ))
                       ds.merge(roots_list[0], other_root)
               # propagate birth elevation to surviving root
               surviving = ds[idx]
               birth_elev[surviving] = min(
                   birth_elev[r] for r in roots_list if not np.isnan(birth_elev[r])
               )
               
       return births, merges, ds
   ```
   
   Call this on the *cropped, subsampled* version of crop_tda 
   (subsample by 4x first: `Z_sub = crop_tda[::4, ::4]`) to keep 
   runtime under 1 minute. Note the resolution change in comments.

3. PRINT STATISTICS:
   - Number of births (local minima found)
   - Number of merge events (saddle points found)
   - Distribution of merge elevations

4. MERGE TREE AS GRAPH: Using networkx, build a tree where:
   - Each node is a basin (identified by its birth elevation + location)
   - Each edge connects two basins at the elevation where they merge
   - Node size = persistence (death - birth)
   
   Draw the tree with elevation on the y-axis (dendrogram style).

Note: this algorithm is O(n log n) for the sort + nearly O(n) for the 
union-find sweep. For the full 3000×3000 grid it may take 30-60 seconds; 
the 4x subsampled version (~750×750) should take ~5 seconds.
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 4 — Priority-queue topological flood simulation

**What it does**: Implements the "hydrologically correct" flood simulation. Starting from river seed pixels, water floods outward in elevation order (using a min-heap). This is connectivity-based — a cell only floods when it is reachable from a seed through already-flooded cells below the current water level. Tracks when separate basins merge.

```
PROMPT 4 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
I want to implement a proper topological flood simulation using a 
priority queue. This differs from the existing `flood_fill` function 
in a key way: instead of BFS from one origin at a fixed threshold, 
this simulates a rising water level starting from all river seed pixels 
simultaneously, expanding only through connected low terrain.

Add a section called "Topological flood simulation (priority queue)".

BACKGROUND:
The algorithm is essentially Prim's / Dijkstra's adapted for terrain 
flooding:
- Start with all seed pixels (from `seed_mask`) added to a min-heap
- The heap stores (elevation, row, col, component_label)
- At each step: pop the lowest-elevation unvisited pixel
- Mark it as flooded with the label of its source component
- Push its unvisited 8-neighbors onto the heap (with their actual elevation)
- When two different-labeled components would merge: record the merge 
  as a saddle event

This is correct because water MUST flow through connected terrain — 
a cell isolated behind a ridge will not flood even if it is low, 
until the water level overtops the ridge.

IMPLEMENTATION:
```python
import heapq

def topological_flood(Z, seed_mask, water_levels=None):
    """
    Z: 2D elevation array (no NaN, barriers as max+1)
    seed_mask: bool array, True = river seed pixels
    water_levels: list of water levels to record snapshots at
                  (default: np.linspace(Z[~(Z>np.nanmax(Z))].min(), 
                                        np.nanpercentile(Z[Z<np.nanmax(Z)], 80), 
                                        50))
    
    Returns:
        flood_labels: 2D array of component labels (-1 = not flooded)
        flood_level_map: 2D array of water level at which each pixel was flooded
        saddle_events: list of (water_level, label_A, label_B)
        snapshots: dict {water_level: flood_labels_copy}
    """
    rows, cols = Z.shape
    flood_labels = np.full((rows, cols), -1, dtype=int)
    flood_level_map = np.full((rows, cols), np.nan)
    
    heap = []
    label_counter = 0
    saddle_events = []
    snapshots = {}
    
    # Initialise from seed pixels
    seed_coords = np.argwhere(seed_mask)
    for r, c in seed_coords:
        label = label_counter
        label_counter += 1
        flood_labels[r, c] = label
        flood_level_map[r, c] = Z[r, c]
        # Push neighbors
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]:
            nr, nc = r+dr, c+dc
            if 0 <= nr < rows and 0 <= nc < cols and flood_labels[nr, nc] < 0:
                heapq.heappush(heap, (Z[nr, nc], nr, nc, label))
    
    # Merge seeds that are already connected (8-neighbors of same initial seeds)
    # Use union-find for label merging
    parent = list(range(label_counter + 10000))  # oversized
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[py] = px
            return True
        return False
    
    if water_levels is None:
        valid = Z[Z < np.nanmax(Z)]
        water_levels = np.linspace(np.nanmin(valid), np.nanpercentile(valid, 80), 50)
    
    water_level_iter = iter(sorted(water_levels))
    next_snapshot = next(water_level_iter, None)
    
    while heap:
        elev, r, c, src_label = heapq.heappop(heap)
        
        if flood_labels[r, c] >= 0:
            # Already flooded — check if different label (saddle event)
            existing = find(flood_labels[r, c])
            incoming = find(src_label)
            if existing != incoming:
                saddle_events.append((elev, existing, incoming))
                union(existing, incoming)
            continue
        
        canonical_label = find(src_label)
        flood_labels[r, c] = canonical_label
        flood_level_map[r, c] = elev
        
        # Snapshot at water levels
        while next_snapshot is not None and elev >= next_snapshot:
            snapshots[next_snapshot] = flood_labels.copy()
            next_snapshot = next(water_level_iter, None)
        
        # Push neighbors
        for dr, dc in [(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(-1,1),(1,-1),(1,1)]:
            nr, nc = r+dr, c+dc
            if 0 <= nr < rows and 0 <= nc < cols and flood_labels[nr, nc] < 0:
                heapq.heappush(heap, (Z[nr, nc], nr, nc, canonical_label))
    
    return flood_labels, flood_level_map, saddle_events, snapshots
```

Run this on `crop_tda` with the `seed_mask` from Module 1.
Use `water_levels = np.linspace(-2, 5, 60)` as a reasonable range 
for the Biesbosch area (NAP-relative, mostly between -2m and +5m).

PRINT:
- Total pixels flooded vs unflooded
- Number of saddle events
- Water level at which each saddle occurs
- Which water level floods 25%, 50%, 75% of the area

PLOT (at the end): 
- `flood_level_map` coloured by the water level at which each cell 
  was flooded (use colormap 'Blues', low levels dark, high levels light)
- Overlay the seed pixels in red
- Title: "Water level at which each cell is first flooded (topological)"
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 5 — Side-by-side comparison (RQ1 + RQ2)

**What it does**: Directly answers RQ1 by comparing the simple threshold model vs. the topological model at the same water levels. Shows isolated depressions that threshold gets wrong.

```
PROMPT 5 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
I now have two flood models:
1. Simple threshold: `simple_flood = crop_tda <= water_level`
2. Topological (from Module 4): `topo_flood = flood_level_map <= water_level`
   (using the flood_level_map from the priority-queue simulation)

Add a section "Comparison: threshold vs topological flood model (RQ1)".

1. SIDE-BY-SIDE SNAPSHOTS: For each of these water levels 
   [-0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0] m NAP:
   Create a 2-column figure (7 rows) showing:
   Left: simple threshold flood extent (blue)
   Right: topological flood extent (blue), seed pixels in red
   Both overlaid on grey terrain. Add water level as the row title.
   Under each pair, print the pixel counts and % difference.

2. OVERESTIMATION MAP: Create a single map showing:
   - Grey: neither model floods this cell
   - Blue: both models agree it's flooded
   - Red: ONLY the threshold model floods it (false positives)
   - Green: ONLY the topological model floods it (shouldn't happen often)
   This directly visualises where the threshold model "cheats" by 
   including isolated low-lying depressions not connected to the river.
   
   Pick a representative water level where the difference is most 
   dramatic (probably around 0.5 – 1.5 m NAP based on the area).

3. FIRST-FLOODED BASINS (RQ2): Using the saddle_events from Module 4,
   identify which basins (by their seed/origin label) first connect to 
   the main river system, and at what water level. 
   Print a ranked table: basin ID, approximate location (centroid), 
   water level at connection.

4. QUANTITATIVE COMPARISON TABLE:
   For water levels [-1, 0, 0.5, 1, 1.5, 2, 3, 5]:
   | Water level (m NAP) | Threshold flooded (%) | Topo flooded (%) | Difference (%) |
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 6 — Merge tree dendrogram + saddle point overlay (RQ3 + RQ4)

**What it does**: Visualizes the merge tree as a dendrogram. Each leaf is a basin; the merge height is the saddle elevation. Directly answers RQ3 (at what elevations do components merge?). Also overlays saddle points on the terrain map.

```
PROMPT 6 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
Using `births`, `merges` from Module 3 (the union-find sweep on the 
subsampled terrain), and `saddle_locations` from Module 2 (GUDHI), 
add a section "Merge tree visualization (RQ3 + RQ4)".

1. DENDROGRAM-STYLE MERGE TREE:
   Build a networkx tree where:
   - Leaf nodes are individual basins (persistent ones, persistence > 0.5m)
   - Internal nodes are merge events
   - Edge weights = difference in elevation between parent and child
   
   Plot this as a dendrogram with:
   - y-axis = elevation (m NAP)
   - Each leaf labeled with its basin's birth elevation and approximate 
     (row, col) location
   - Internal nodes labeled with their merge elevation
   - Nodes colored by persistence (use a colormap, high persistence = vivid)
   - Horizontal dashed lines at key water levels (e.g., 0m, 0.5m, 1m, 2m)
     showing which components have merged by that level
   
   This is the "topological fingerprint" of the terrain.

2. SADDLE POINT MAP: On the terrain image (crop_tda), overlay:
   - Basin minima: blue circles, sized by persistence
   - Saddle points: red triangles (from GUDHI saddle_locations)
   - River seeds: cyan dots
   - Draw thin lines connecting each minimum to its saddle
   This directly answers RQ4: these structures (barriers between basins) 
   are invisible from raw elevation.

3. PERSISTENCE-FILTERED COMPARISON: Show three versions of the basin 
   minima map side by side:
   - All local minima (very noisy)
   - Filtered by persistence > 0.1m (still noisy)
   - Filtered by persistence > 0.5m (clean, meaningful basins)
   Caption: "The persistence threshold determines which depressions are 
   real vs. LiDAR noise"

4. CRITICAL ELEVATION TABLE: Print a table of all merge events above 
   persistence threshold 0.5m, sorted by merge elevation:
   | Rank | Merge elevation (m NAP) | Basin A birth elev | Basin B birth elev | Persistence |
   These are the "critical thresholds" of RQ3.
───────────────────────────────────────────────────────────────────────────────
```

---

### MODULE 7 — Flood animation

**What it does**: Generates an animated GIF / interactive animation of the topological flood simulation. Colored by connected component (which basin). This is the main visual deliverable.

```
PROMPT 7 (paste into Claude Code)
───────────────────────────────────────────────────────────────────────────────
Using the `snapshots` dict from Module 4 (keys = water levels, 
values = flood_labels arrays), create an animation of the topological 
flood simulation.

Add a section "Flood animation".

1. COMPONENT COLORING: Before animating, assign a consistent color to 
   each component label. Use a qualitative colormap (e.g., 'tab20') 
   with enough colors for the expected number of components. 
   Unflooded pixels: white/transparent. 
   Flooded pixels: colored by their component label.
   Note: component labels change as merges happen (union-find merges 
   them into one label). Use the canonical label from the parent array.

2. ANIMATION:
   Use `matplotlib.animation.FuncAnimation`. Each frame:
   - Background: terrain in grey ('Greys_r', low alpha)
   - Flooded pixels: colored by component, full alpha
   - Title: f"Water level: {water_level:.2f} m NAP | Flooded: {pct:.1f}%"
   - Red markers for seed pixels (constant across frames)
   
   Frames: use all keys in `snapshots` sorted by water level.
   Interval: 150ms between frames.
   
   Save as both:
   - `flood_animation.gif` (for the report, ~60 frames)
   - `flood_animation.mp4` (higher quality, if ffmpeg available)

3. KEY-FRAME FIGURE: Extract 6 frames at evenly spaced water levels 
   and plot them as a 2×3 figure. This is the static version for the 
   report.

4. SIMPLE THRESHOLD PARALLEL: Make a second animation with identical 
   framing but showing the simple threshold model (all cells below 
   water level, no connectivity, single blue color). 
   Save as `flood_threshold_animation.gif`.
   
   The visual contrast between the two animations is the clearest 
   illustration of RQ1.

Note: for performance, downsample both `crop_tda` and `flood_level_map` 
by 2x before animating if the 3000×3000 full resolution is too slow.
───────────────────────────────────────────────────────────────────────────────
```

---

## 5. Execution order

Run the modules in order 1 → 7. Each module depends on variables from previous ones. Suggested checkpoints:

```
Module 1 → produces: seed_mask, crop_tda
Module 2 → produces: persistence_pairs, basin_locations, saddle_locations
Module 3 → produces: births, merges (on subsampled grid)
Module 4 → produces: flood_labels, flood_level_map, saddle_events, snapshots
Module 5 → produces: comparison figures, overestimation map
Module 6 → produces: merge tree dendrogram, saddle point map
Module 7 → produces: flood_animation.gif, key-frame figure
```

---

## 6. Connecting back to the report (TDA techniques section)

For each figure you produce, the report section maps to a research question:

| Figure | Report section | RQ answered |
|---|---|---|
| Persistence diagram (Module 2) | TDA Techniques → GUDHI, sublevel set persistence | Motivation for method |
| Overestimation map (Module 5) | Results → Threshold vs. topological | RQ1 |
| First-flooded basins table (Module 5) | Results → Basin connectivity | RQ2 |
| Merge tree dendrogram + critical elevations (Module 6) | Results → Saddle points | RQ3 |
| Saddle point map (Module 6) | Results → Topological structures | RQ4 |
| Animation key frames (Module 7) | Discussion / Appendix | All RQs visually |

---

## 7. Common pitfalls to flag for Claude Code

When using these prompts, watch for:

- **Memory**: The full 3000×3000 grid has 9M pixels. GUDHI CubicalComplex on this may use ~2GB RAM. Always subsample 2-4x for the merge tree sweep; use the full grid for visualization only.
- **GUDHI's NaN behavior**: Pass `crop_tda` (no NaN) not `crop`. GUDHI does not handle NaN gracefully.
- **DisjointSet indexing**: `scipy.cluster.hierarchy.DisjointSet` uses objects as keys, not integers. Use pixel indices (integers) and wrap in `DisjointSet(range(N))`.
- **Union-find in the priority queue (Module 4)**: After a merge, all pixels that were labeled `A` still have label `A` in the array. Use `find(label)` everywhere instead of the raw label to get the canonical component.
- **Animation color consistency**: Component label 0 might be the biggest component (the main river) — make sure it gets a stable, intuitive color (e.g., light blue) across all frames.
