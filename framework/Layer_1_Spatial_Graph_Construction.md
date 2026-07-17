
# MOSAIC Layer 1 — Spatial Graph Construction

Biological Entities

↓

Spatial Graph



## Purpose

Layer 1 transforms the biological entities defined in Layer 0 into a spatial graph.

Each detected cell becomes a node, and spatial relationships between cells become edges. This transformation converts a cell-by-cell coordinate table into a mathematical representation of tissue architecture.

The central question of Layer 1 is:

> How does tissue become a graph?

The resulting graph provides the structural foundation for all downstream MOSAIC analyses.

---

## Position Within the MOSAIC Framework

```text
Layer 0: Biological Entities
            ↓
Layer 1: Spatial Graph Construction
            ↓
Layer 2: Network Features
```

Layer 0 defines the biological objects entering the framework.

Layer 1 defines how those objects are connected.

Layer 2 measures the resulting network.

---

## Inputs from Layer 0

Layer 1 receives a standardized cell-level dataset from Layer 0.

At minimum, each cell should include:

* Cell ID
* X coordinate
* Y coordinate
* Cell phenotype or cell type
* Sample or slide identifier

Optional metadata may include:

* Tissue compartment
* Disease type
* Patient identifier
* Clinical outcome
* Imaging platform
* Marker intensities
* Gene expression
* Protein abundance
* Mutation status
* Cell quality metrics
* Time point
* Experimental batch

The graph construction process should preserve the original cell identifiers and metadata so that network results can be traced back to the biological source data.

---

## Nodes

Each biological entity becomes a graph node.

In the standard MOSAIC implementation, the primary node is the individual cell.

A cell node may contain attributes such as:

```text
Cell ID
X coordinate
Y coordinate
Cell phenotype
Sample ID
Tissue compartment
Molecular measurements
Clinical metadata
```

The node represents the biological entity itself.

Its attributes describe what the entity is, while its edges describe how it is positioned relative to other entities.

Although cells are the default graph unit, MOSAIC can also support other biological entities, including:

* Tissue regions
* Spatial spots
* Cell neighborhoods
* Glands
* Vessels
* Tumor nests
* Tertiary lymphoid structures
* Image-derived objects
* Molecular features

This allows the framework to support both single-scale and multiscale tissue representations.

---

## Edges

Edges represent spatial relationships between biological entities.

In the simplest implementation, two cells are connected when they fall within a predefined spatial neighborhood.

An edge may represent:

* Physical proximity
* Direct adjacency
* Shared neighborhood membership
* Potential cell-cell interaction
* Anatomical continuity
* Similarity in molecular or phenotypic state
* A user-defined biological relationship

For a radius-based graph, an edge is created when:

```text
distance(cell i, cell j) ≤ radius threshold
```

The edge indicates that the two cells occupy the same local spatial environment under the selected neighborhood definition.

An edge does not necessarily prove that two cells are communicating biologically. It records a measurable spatial relationship that can later be evaluated in biological context.

---

## Neighborhood Definitions

The definition of a neighborhood is one of the most important decisions in spatial graph construction.

MOSAIC does not require a single universal neighborhood method. The appropriate strategy depends on the imaging platform, tissue scale, cell density, biological question, and spatial resolution.

Supported approaches may include the following.

---

### Radius-Based Graph

In a radius-based graph, two cells are connected when the Euclidean distance between their coordinates is less than or equal to a specified threshold.

```text
Cell A ─── Cell B

Connected if distance ≤ R
```

A radius graph is biologically intuitive because the threshold can be selected to approximate a meaningful interaction distance.

Advantages include:

* Direct interpretation in physical units
* Efficient construction using spatial indexing methods such as KDTree
* Natural representation of local cellular neighborhoods
* Compatibility with tissues containing irregular cell distributions

Limitations include:

* Sensitivity to the selected radius
* Variable node degree in areas with different cell densities
* Potential differences across imaging platforms or tissue types

MOSAIC currently uses radius-based graph construction as a primary reference implementation.

A threshold such as 20 µm may be appropriate for direct local cellular proximity in some datasets, but the optimal value should be evaluated for each biological and technical context.

---

### k-Nearest Neighbors

In a k-nearest-neighbor graph, each cell is connected to a fixed number of its closest neighboring cells.

```text
Each cell → k closest cells
```

Advantages include:

* Every cell receives a defined number of neighbors
* Reduced sensitivity to local density variation
* Useful for datasets in which physical scale is uncertain or inconsistent

Limitations include:

* Cells may be connected across biologically large distances in sparse regions
* The same value of k may represent different physical distances across a sample
* Directed and undirected implementations may produce different graph structures

A k-nearest-neighbor graph may be appropriate when relative neighborhood structure is more important than an absolute physical interaction radius.

---

### Delaunay Triangulation

Delaunay triangulation constructs edges by connecting spatial points according to their geometric arrangement.

This approach creates a graph in which neighboring points are connected without requiring a fixed radius or fixed number of neighbors.

Advantages include:

* Geometry-driven neighborhood construction
* Adaptation to local cell density
* Natural representation of spatial adjacency

Limitations include:

* Possible long edges near tissue boundaries
* Sensitivity to sparse regions and segmentation artifacts
* Less direct biological interpretation than a fixed physical radius

Delaunay graphs can be useful when tissue geometry and local packing are central to the analysis.

---

### Voronoi-Based Adjacency

Voronoi tessellation divides space into regions surrounding each cell coordinate.

Two cells may be considered neighbors when their Voronoi regions share a boundary.

This approach approximates direct spatial adjacency and can be useful when cell packing and territorial relationships are important.

However, Voronoi-based approaches may be affected by tissue boundaries, missing cells, segmentation quality, and highly irregular cell distributions.

---

### Custom Biological Adjacency

Some biological questions require a relationship definition that extends beyond geometric proximity.

Custom edges may represent:

* Tumor cells connected to nearby immune cells
* Cells belonging to the same gland
* Cells sharing a vascular structure
* Cells within the same tertiary lymphoid structure
* Ligand-receptor-compatible cell pairs
* Cells linked through a defined anatomical compartment
* Cells observed at consecutive time points

Custom adjacency allows the graph structure to reflect the biological hypothesis being tested.

Such edges should remain clearly distinguished from purely spatial edges.

---

## Distance Calculation

For two cells with coordinates:

```text
Cell i = (xi, yi)
Cell j = (xj, yj)
```

the Euclidean distance is:

```text
d(i,j) = √[(xi - xj)² + (yi - yj)²]
```

For three-dimensional data, the z coordinate can be included:

```text
d(i,j) = √[(xi - xj)² + (yi - yj)² + (zi - zj)²]
```

Coordinate units must be documented clearly.

Examples include:

* Micrometers
* Pixels
* Normalized coordinates
* Spatial transcriptomics spot positions

Pixel coordinates should be converted to physical units when reliable image calibration information is available.

---

## Graph Types

The graph representation may vary according to the biological question.

---

### Undirected Graph

In an undirected graph, the relationship between two cells is symmetric.

```text
Cell A ─── Cell B
```

If Cell A is near Cell B, then Cell B is also near Cell A.

Undirected graphs are appropriate for most physical proximity analyses.

---

### Directed Graph

In a directed graph, the relationship has a defined orientation.

```text
Cell A ───> Cell B
```

Directed edges may be useful for:

* Predicted signaling direction
* Temporal transitions
* Cell migration
* Differentiation trajectories
* Source-target biological relationships

Direction should not be inferred from spatial proximity alone.

---

### Weighted Graph

Edges may include numerical weights.

Possible edge weights include:

* Euclidean distance
* Inverse distance
* Interaction probability
* Ligand-receptor score
* Contact strength
* Shared-boundary length
* Confidence score

For example:

```text
weight(i,j) = 1 / distance(i,j)
```

In this formulation, closer cells receive stronger edge weights.

The weighting method should be selected according to the biological interpretation intended for downstream analysis.

---

### Multilayer or Multiplex Graph

A multilayer graph may contain several edge types between the same biological entities.

For example:

```text
Cell A ── spatial proximity ── Cell B
Cell A ── ligand-receptor ──── Cell B
Cell A ── shared community ─── Cell B
```

This allows physical, molecular, anatomical, and functional relationships to coexist within the same graph framework.

Multilayer graphs are particularly relevant to the later multiomic stages of MOSAIC.

---

### Bipartite or Heterogeneous Graph

MOSAIC may also represent different categories of nodes.

Examples include:

```text
Cell ── ANNOTATED_AS ──> Cell Type

Cell ── BELONGS_TO ────> Sample

Cell ── IN_CLUSTER ────> Community

Cell ── LOCATED_IN ────> Tissue Region
```

This structure is especially useful for knowledge-graph implementations such as Neo4j.

---

## Reference Graph Schema

A basic MOSAIC graph may contain the following structure:

```text
(:Cell)-[:NEAR_TO]-(:Cell)
```

Each `Cell` node may contain properties such as:

```text
cell_id
x
y
phenotype
sample_id
compartment
```

Each `NEAR_TO` relationship may contain properties such as:

```text
distance
radius
weight
construction_method
```

A heterogeneous implementation may also include:

```text
(:Cell)-[:ANNOTATED_AS]->(:CellType)

(:Cell)-[:BELONGS_TO_SAMPLE]->(:Sample)

(:Cell)-[:IN_CLUSTER]->(:Community)
```

This schema preserves both the cell-level spatial network and the biological context required for interpretation.

---

## Sample Boundaries

Graphs should normally be constructed independently for each biological sample, slide, tissue section, or region of interest.

Cells from unrelated samples must not be connected solely because their coordinate values happen to overlap.

A typical workflow is:

```text
Split dataset by sample
        ↓
Construct one graph per sample
        ↓
Calculate sample-level features
        ↓
Compare graphs across samples
```

This preserves the independence of biological observations and prevents artificial cross-sample edges.

---

## Tissue Boundaries and Compartments

Spatial graph construction should consider tissue boundaries and anatomical compartments.

Potential issues include:

* Connecting cells across empty space
* Connecting cells across tissue folds
* Connecting cells across separate tissue fragments
* Connecting cells across luminal spaces
* Connecting cells across manually defined regions
* Connecting cells across tumor-stroma boundaries when such edges are not biologically intended

When available, masks, tissue regions, segmentation boundaries, or compartment labels may be used to constrain graph construction.

For example, an edge may require both:

```text
distance ≤ R
```

and

```text
same tissue section or permitted compartment relationship
```

This prevents geometrically close but biologically separated cells from being treated as direct neighbors.

---

## Quality Control

Graph quality depends directly on the quality of the Layer 0 input data.

Before graph construction, the following should be evaluated:

* Missing coordinates
* Duplicate cell identifiers
* Duplicate coordinates
* Invalid sample identifiers
* Segmentation artifacts
* Outlier coordinates
* Incorrect coordinate units
* Tissue fragments
* Low-quality cell detections
* Extreme differences in cell density
* Phenotype annotation quality

After graph construction, useful quality-control summaries include:

* Number of nodes
* Number of edges
* Mean and median degree
* Isolated-node count
* Connected-component count
* Edge-distance distribution
* Cell-density distribution
* Graph density
* Sample-specific graph size

Unexpected values may indicate an unsuitable neighborhood threshold, coordinate error, or segmentation problem.

---

## Radius Sensitivity

The selected spatial radius influences the resulting graph topology.

A small radius may produce:

* Many isolated cells
* Fragmented components
* Highly local interactions
* Reduced connectivity

A large radius may produce:

* Dense graphs
* Loss of local specificity
* Increased computational cost
* Connections between cells unlikely to interact directly

MOSAIC therefore encourages sensitivity analysis across biologically plausible radii.

For example:

```text
10 µm
20 µm
30 µm
40 µm
```

The goal is not necessarily to identify one universally correct radius, but to determine whether biological conclusions remain stable across reasonable spatial scales.

Multiscale graphs may also be constructed when tissue organization is expected to operate at several distances.

---

## Computational Implementation

Efficient graph construction is essential for large spatial datasets.

A naive all-versus-all distance comparison has approximately quadratic computational complexity and becomes impractical for large cell populations.

Spatial indexing methods such as KDTree can identify nearby cells efficiently.

A conceptual workflow is:

```text
Load cell coordinates
        ↓
Separate cells by sample
        ↓
Build spatial index
        ↓
Query neighboring cells
        ↓
Create node and edge tables
        ↓
Export graph
```

Possible output formats include:

* NetworkX graph objects
* igraph objects
* Edge-list CSV files
* Node-list CSV files
* GraphML
* GEXF
* Neo4j import tables
* AnnData-compatible structures

The graph representation should remain interoperable so that users can analyze it using different graph and network-science platforms.

---

## Example Workflow

Assume a cell-level table contains:

```text
cell_id,x,y,cell_type,sample_id
C001,120.4,84.2,T_cell,S01
C002,132.1,91.6,Fibroblast,S01
C003,205.5,140.7,Tumor,S01
C004,127.8,79.9,Macrophage,S01
```

Using a radius threshold of 20 µm:

1. Each row becomes a `Cell` node.
2. Pairwise spatial distances are evaluated using a spatial index.
3. Cells separated by 20 µm or less are connected.
4. The edge distance is stored as a relationship property.
5. The resulting graph is exported for Layer 2 analysis.

A simplified graph may appear as:

```text
T Cell ───── Fibroblast
   │             │
   └──── Macrophage

Tumor Cell
```

The isolated tumor cell is still retained as a node, even if it has no neighbors within the selected radius.

This absence of local connections may itself carry biological information.

---

## Outputs

Layer 1 produces a spatial graph for each sample.

The primary outputs are:

### Node Table

```text
node_id
cell_id
x
y
cell_type
sample_id
additional metadata
```

### Edge Table

```text
source
target
distance
weight
edge_type
sample_id
```

### Graph Object

A computational graph that can be analyzed using graph-theory, network-science, machine-learning, or knowledge-graph tools.

### Quality-Control Summary

```text
number of nodes
number of edges
mean degree
isolated cells
connected components
edge-distance distribution
```

These outputs become the input for Layer 2.

---

## Design Principles

### Spatial relationships must be explicit

MOSAIC requires the rule used to create an edge to be documented clearly.

A graph should never be constructed using an undefined or hidden neighborhood assumption.

### Construction is separate from interpretation

Layer 1 records relationships.

It does not claim that proximity proves communication, causality, or functional interaction.

### Biological scale must guide graph scale

The neighborhood definition should reflect the spatial scale relevant to the biological question.

### Sample identity must be preserved

Graphs should not connect independent samples unless a separate cross-sample analytical graph is intentionally constructed.

### Metadata must remain traceable

Every node and edge should be traceable to the original dataset and graph-construction parameters.

### Multiple graph definitions may be valid

Different neighborhood strategies may capture different aspects of tissue architecture.

MOSAIC treats graph construction as a configurable analytical decision rather than a fixed universal rule.

### Reproducibility is essential

The following should be recorded:

* Neighborhood method
* Radius or k value
* Coordinate units
* Edge-weighting method
* Tissue-boundary constraints
* Software version
* Graph-construction date
* Input dataset version

---

## Relationship to Layer 2

Layer 1 creates the graph.

Layer 2 measures it.

```text
Biological Entities
        ↓
Spatial Graph
        ↓
Degree
Clustering
Centrality
Communities
Topology
```

The same cells may generate very different network structures depending on how they are spatially organized.

This is the central transition from cell composition to tissue architecture.

---

## Summary

Layer 1 converts spatially resolved biological entities into a graph representation.

Cells become nodes.

Spatial relationships become edges.

The resulting graph preserves the architecture of the tissue while remaining independent of downstream biological interpretation.

Layer 1 therefore establishes the computational substrate on which MOSAIC measures network features, identifies organizational states, and constructs tissue fingerprints.

> Layer 1 builds biological topology. It does not yet assign biological meaning.
