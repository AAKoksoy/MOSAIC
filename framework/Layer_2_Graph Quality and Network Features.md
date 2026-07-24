

# MOSAIC Layer 2 — Graph Quality and Network Features

## Purpose

Layer 2 evaluates whether a spatial graph is technically credible and measures the organizational properties that emerge from it.

Layer 1 converts biological observations into a graph by representing cells as nodes and spatial relationships as edges.

Layer 2 asks two related questions:

1. **Is the constructed graph trustworthy?**
2. **What aspects of tissue organization can be measured from that graph?**

For this reason, Layer 2 is divided into two analytical stages:

* **Layer 2A — Graph Quality Assessment**
* **Layer 2B — Biological Network Features**

Layer 2A

Graph Quality Assessment

↓

Layer 2B

Biological Network Features

Graph metrics should not be interpreted biologically until the graph itself has been evaluated for technical artifacts, inappropriate assumptions, and construction errors.

---

## Scientific Question

**Is the graph trustworthy, and what measurable organizational properties emerge from it?**

---

## Inputs from Layer 1

Layer 2 receives the outputs of Layer 1.

The minimum required inputs are:

### Cell table

Each row represents one biological entity, initially a cell.

Required fields include:

```text
cell_id
x
y
cell_type
sample_id
```

Optional fields may include:

```text
patient_id
disease_type
tissue_compartment
marker_expression
gene_expression
clinical_outcome
quality_control_metrics
```

### Spatial edge table

Each row represents one spatial relationship.

Recommended fields include:

```text
source
target
distance_um
sample_id
relationship_type
```

The initial spatial relationship type is:

```text
NEAR_TO
```

The `NEAR_TO` relationship indicates that two cells satisfy a defined spatial-neighborhood rule.

It represents a spatial proposition, not proof of a biological interaction.

---

# Layer 2A — Graph Quality Assessment

## Why Graph Quality Must Be Evaluated

A graph is not automatically biologically meaningful simply because it was successfully constructed.

Graph topology can be altered by:

* incorrect coordinate units;
* inappropriate neighborhood thresholds;
* segmentation artifacts;
* duplicated cells;
* missing tissue regions;
* slide-border effects;
* cross-sample connections;
* unusually dense cell regions;
* image-registration errors;
* inconsistent cell detection;
* self-loops or duplicate edges.

These issues can produce network patterns that appear biologically interesting but are actually technical artifacts.

MOSAIC therefore treats graph quality assessment as a required analytical stage rather than an optional preprocessing step.

---

## Core Graph Quality Metrics

### Node Count

The number of cells represented in the graph.

```text
number_of_nodes
```

The graph node count should generally agree with the number of validated cells received from Layer 0.

A difference between the validated cell count and graph node count may indicate:

* cells removed during graph construction;
* isolated cells omitted from the graph object;
* duplicate identifiers;
* filtering that was not recorded;
* incomplete graph creation.

---

### Edge Count

The number of spatial relationships present in the graph.

```text
number_of_edges
```

Edge count depends strongly on:

* the selected neighborhood method;
* the radius or nearest-neighbor parameter;
* cell density;
* tissue composition;
* sample size;
* image resolution.

Edge count should therefore be interpreted in the context of the graph-construction assumptions defined in Layer 1.

---

### Graph Density

Graph density describes the proportion of all possible node pairs that are connected.

For an undirected graph:

```text
density = 2E / N(N - 1)
```

where:

```text
E = number of edges
N = number of nodes
```

Spatial graphs are usually sparse because cells connect only to local neighbors.

An unexpectedly dense graph may indicate:

* an excessively large radius;
* incorrect coordinate scaling;
* duplicated spatial positions;
* unintended cross-sample edges.

An unexpectedly sparse graph may indicate:

* an excessively small radius;
* incorrect spatial units;
* missing coordinates;
* sparse tissue sampling;
* over-filtering.

Graph density should not be compared across samples without considering differences in tissue area, cell abundance, and acquisition platform.

---

### Degree Distribution

The degree of a node is the number of edges connected to that node.

For a cell node:

```text
degree(cell) = number of spatial neighbors
```

The graph-wide degree distribution provides an overview of local connectivity.

Important summaries include:

```text
minimum degree
maximum degree
mean degree
median degree
degree standard deviation
degree quantiles
```

Degree distribution can reveal:

* isolated cells;
* highly crowded regions;
* extreme connectivity outliers;
* possible segmentation artifacts;
* differences in cell density;
* inappropriate graph parameters.

Degree is therefore both a quality-control measure and a biological network feature.

---

### Isolated Nodes

An isolated node has degree zero.

```text
degree = 0
```

Possible explanations include:

* a genuinely isolated cell;
* a cell outside the selected radius of all other cells;
* a tissue-border cell;
* debris or false-positive segmentation;
* incorrect coordinates;
* an overly restrictive neighborhood threshold.

Isolated nodes should not automatically be removed.

They should first be reported and evaluated because isolation may represent either a technical artifact or a biologically meaningful state.

---

### High-Degree Outliers

A high-degree outlier is a node with substantially more neighbors than most other nodes in the graph.

Possible explanations include:

* a genuinely dense cellular region;
* a highly connected biological hub;
* an excessively large segmented object;
* duplicated cells;
* coordinate overlap;
* incorrect coordinate units;
* a graph-construction error.

High-degree cells should be identified using transparent and configurable criteria.

Possible approaches include:

```text
interquartile-range threshold
percentile threshold
z-score threshold
median absolute deviation
```

No single outlier threshold should be treated as universally correct.

---

### Self-Loops

A self-loop occurs when a node is connected to itself.

```text
source == target
```

Self-loops should normally be excluded from the initial undirected spatial graph.

Their presence may indicate an error during:

* nearest-neighbor querying;
* pair generation;
* edge-table construction;
* graph conversion.

---

### Duplicate Edges

A duplicate edge occurs when the same cell pair is represented more than once.

For an undirected graph:

```text
Cell_A — Cell_B
```

and:

```text
Cell_B — Cell_A
```

represent the same relationship unless the graph is intentionally directed.

Duplicate edges can artificially inflate:

* degree;
* weighted degree;
* graph density;
* interaction counts;
* centrality scores.

Layer 2 should verify that each undirected cell pair appears only once.

---

### Cross-Sample Edges

Cells from different samples should not normally be connected through spatial proximity.

A cross-sample edge occurs when:

```text
source.sample_id != target.sample_id
```

Possible causes include:

* graph construction across a combined dataset without sample grouping;
* duplicated coordinates across slides;
* missing sample identifiers;
* incorrect data merging.

Cross-sample edges should be treated as graph-construction errors unless they were explicitly created for a separate analytical purpose.

---

### Edge-Distance Distribution

Each spatial edge should retain the Euclidean distance between its source and target cells.

Important summaries include:

```text
minimum distance
maximum distance
mean distance
median distance
distance quantiles
```

For a radius graph, no edge distance should exceed the selected radius, allowing for minor numerical tolerance.

Unexpected distance values may reveal:

* incorrect coordinate units;
* errors in distance calculation;
* unintended graph construction;
* mixed pixel and micrometer measurements;
* incorrect metadata.

---

### Connected Components

A connected component is a set of nodes in which every node can be reached from every other node through a sequence of edges.

Graph-quality summaries should include:

```text
number of connected components
size of largest component
size of smallest component
median component size
fraction of nodes in largest component
number of singleton components
```

A graph containing one giant component may represent continuous tissue organization, but it may also result from an overly permissive neighborhood threshold.

A graph containing many small components may represent fragmented tissue, but it may also indicate an overly restrictive threshold.

Connected-component structure must therefore be evaluated in relation to the tissue, acquisition platform, and scientific question.

---

### Tissue-Border Effects

Cells close to an image boundary may have fewer observable neighbors because surrounding tissue falls outside the imaged region.

This can lower:

* degree;
* clustering coefficient;
* local density;
* neighborhood diversity.

Low connectivity at a tissue or image border should not automatically be interpreted as biological isolation.

Future MOSAIC implementations may incorporate:

* distance to image boundary;
* tissue-mask boundaries;
* field-of-view boundaries;
* border-cell flags;
* edge-corrected neighborhood statistics.

---

### Segmentation Artifacts

Graph topology can expose possible segmentation problems.

Examples include:

* extremely high-degree objects representing merged cells;
* groups of cells with identical coordinates;
* isolated objects representing debris;
* artificial gaps caused by under-segmentation;
* unusually dense clusters caused by duplicated detections.

Layer 2 does not replace image-level quality control.

Instead, it provides a graph-based view of potential segmentation abnormalities that may not be obvious from the cell table alone.

---

## Recommended Graph Quality Report

A minimum graph-quality report should include:

```text
number_of_nodes
number_of_edges
graph_density
number_of_isolated_nodes
fraction_isolated
number_of_self_loops
number_of_duplicate_edges
number_of_cross_sample_edges
number_of_connected_components
largest_component_size
fraction_in_largest_component
minimum_degree
maximum_degree
mean_degree
median_degree
minimum_edge_distance
maximum_edge_distance
mean_edge_distance
median_edge_distance
```

Graph-quality results should be generated for:

* the complete dataset;
* each sample separately;
* optionally each tissue compartment.

---

# Layer 2B — Biological Network Features

## Purpose

After graph quality has been assessed, network features can be used to describe tissue organization.

These metrics measure different aspects of the structural position of each cell.

The initial MOSAIC implementation includes:

* degree;
* weighted degree;
* PageRank;
* betweenness centrality;
* clustering coefficient;
* connected-component membership.

Each metric should be interpreted as a structural measurement rather than direct proof of biological function.

---

## Degree

### Definition

Degree is the number of graph edges connected to a node.

For a cell in an undirected spatial graph:

```text
degree(cell) = number of neighboring cells
```

### Unit of Analysis

Degree is calculated per cell.

It can later be summarized by:

* cell type;
* tissue compartment;
* sample;
* patient;
* disease group;
* clinical outcome.

### Possible Biological Interpretation

Degree may reflect:

* local cellular crowding;
* neighborhood exposure;
* spatial integration;
* possible hub status;
* the number of opportunities for local interaction.

### Limitations

Degree depends directly on:

* the selected spatial radius;
* the neighborhood method;
* cell density;
* segmentation quality;
* tissue boundaries;
* coordinate scaling.

A high degree does not necessarily indicate biological importance.

Abundant cell types may have high average degree simply because they occupy dense tissue regions.

---

## Weighted Degree

### Definition

Weighted degree is the sum of the weights of all edges connected to a node.

```text
weighted_degree(i) = Σ weight(i, j)
```

The biological meaning depends entirely on how edge weights are defined.

Possible edge weights include:

* inverse distance;
* normalized proximity;
* interaction confidence;
* ligand–receptor evidence;
* marker-expression similarity;
* repeated-observation strength.

### Distance-Based Weighting

For spatial graphs, a possible weight is:

```text
weight = 1 / (distance + ε)
```

where `ε` prevents division by zero.

Under this definition, closer neighbors contribute more strongly than distant neighbors.

### Possible Biological Interpretation

Weighted degree may reflect:

* total local relationship strength;
* cumulative proximity;
* intensity of neighborhood embedding;
* strength of evidence supporting local connections.

### Limitations

Weighted degree cannot be interpreted without documenting the edge-weight definition.

Different weighting functions can produce different conclusions from the same graph.

MOSAIC must preserve both:

```text
distance_um
edge_weight
```

when weighted metrics are used.

---

## PageRank

### Definition

PageRank estimates the structural importance of a node by considering both:

* how many connections the node has;
* how structurally important its neighbors are.

A cell connected to several highly connected cells may receive a higher PageRank score than a cell with the same degree connected mainly to peripheral cells.

### Unit of Analysis

PageRank is calculated per cell.

### Possible Biological Interpretation

PageRank may reflect:

* network influence;
* integration into highly connected neighborhoods;
* structural importance beyond raw neighbor count;
* participation in dominant tissue organization.

### Limitations

PageRank was not originally designed to measure biological influence.

Its interpretation depends on:

* graph directionality;
* edge weighting;
* damping factor;
* graph density;
* component structure.

A high PageRank score does not establish signaling activity, regulatory control, or causal biological influence.

It identifies structural prominence within the defined graph.

---

## Betweenness Centrality

### Definition

Betweenness centrality measures how often a node lies along shortest paths between other nodes.

Conceptually:

```text
betweenness(v) =
fraction of shortest paths between node pairs that pass through v
```

### Unit of Analysis

Betweenness is calculated per cell.

### Possible Biological Interpretation

High-betweenness cells may represent:

* bridges between spatial regions;
* interface cells;
* transition zones;
* connectors between cellular communities;
* possible communication bottlenecks.

### Limitations

Betweenness can be sensitive to:

* graph size;
* disconnected components;
* local graph density;
* radius selection;
* small changes in graph construction.

Exact betweenness calculations may also be computationally expensive for large spatial graphs.

Approximation methods may be required for large datasets and must be clearly documented.

A high-betweenness cell is a structural bridge in the graph. It is not automatically a biological mediator.

---

## Clustering Coefficient

### Definition

The local clustering coefficient measures how frequently a node’s neighbors are connected to one another.

For an undirected graph:

```text
clustering_coefficient =
observed connections among neighbors
/
possible connections among neighbors
```

### Unit of Analysis

Clustering coefficient is calculated per cell.

### Possible Biological Interpretation

A high clustering coefficient may indicate:

* compact local organization;
* cohesive cellular neighborhoods;
* repeated local contacts;
* structured microenvironments;
* membership in a tightly organized tissue region.

A low clustering coefficient may indicate:

* a peripheral cell;
* a bridge between regions;
* an interface position;
* a sparse or linear arrangement.

### Limitations

Clustering coefficient depends on:

* degree;
* graph-construction method;
* radius;
* tissue density;
* tissue boundaries.

Cells with fewer than two neighbors generally have a clustering coefficient of zero because no neighbor-neighbor triangle can be formed.

Low values in low-degree or border cells should therefore be interpreted cautiously.

---

## Connected-Component Membership

### Definition

Connected-component analysis assigns each node to an uninterrupted graph structure.

Each cell receives:

```text
component_id
component_size
```

### Possible Biological Interpretation

Connected components may identify:

* isolated cells;
* small cellular aggregates;
* fragmented tissue regions;
* large continuous tissue systems;
* spatial structures created by hypothesis-specific relationships.

### Limitations

In a raw proximity graph, connected components are determined largely by the selected neighborhood threshold.

They should not automatically be interpreted as biological communities.

Connected components become especially informative when applied to a hypothesis-specific graph containing biologically qualified relationships.

---

# Cell-Level Output Table

The initial Layer 2 implementation should generate one row per cell.

Recommended fields include:

```text
cell_id
sample_id
cell_type
degree
weighted_degree
pagerank
betweenness
clustering_coefficient
component_id
component_size
```

Optional fields may include:

```text
is_isolated
is_high_degree_outlier
distance_to_boundary
tissue_compartment
patient_id
disease_type
```

The original biological metadata should remain linked to the calculated network features.

---

# Aggregation by Cell Type

Cell-level network metrics may be summarized by phenotype or cell type.

Recommended summaries include:

```text
cell_count
mean
median
standard_deviation
minimum
maximum
25th_percentile
75th_percentile
```

Example output:

```text
sample_id
cell_type
cell_count
mean_degree
median_degree
mean_pagerank
median_pagerank
mean_betweenness
mean_clustering_coefficient
```

Cell-type aggregation can reveal whether particular cell populations tend to occupy:

* highly connected regions;
* cohesive neighborhoods;
* bridge positions;
* peripheral positions;
* structurally influential regions.

However, cell abundance must be considered during interpretation.

A common cell type may dominate raw network summaries because it contributes more nodes and edges.

MOSAIC should distinguish between:

* raw network contribution;
* per-cell average;
* abundance-normalized contribution;
* sample-level prevalence.

---

# Aggregation by Sample

Network features should also be summarized independently for each sample.

Recommended sample-level summaries include:

```text
sample_id
number_of_cells
number_of_edges
graph_density
mean_degree
median_degree
mean_clustering_coefficient
number_of_components
largest_component_fraction
mean_pagerank
mean_betweenness
```

Sample-level summaries provide a basis for:

* patient comparison;
* disease-subtype comparison;
* treatment-response analysis;
* tissue-fingerprint construction;
* quality-control review.

Metrics should not be compared across samples without considering differences in:

* cell count;
* tissue area;
* sampling density;
* platform;
* segmentation;
* neighborhood parameters;
* tissue compartment composition.

---

# Raw Metrics and Normalized Metrics

MOSAIC should preserve raw network metrics while also allowing normalized comparisons.

Possible normalization strategies include:

* division by sample mean;
* z-score normalization within sample;
* percentile rank within sample;
* abundance normalization by cell type;
* comparison against label permutations;
* comparison against spatial null models.

Normalization should never replace the raw metric.

Both the original measurement and the transformation should be retained.

Example:

```text
degree
degree_zscore_within_sample
degree_percentile_within_sample
degree_abundance_normalized
```

The selected normalization strategy must be documented because different strategies answer different scientific questions.

---

# Graph Metrics Are Not Biological Conclusions

Network metrics describe the topology of the graph that was constructed.

They do not independently prove:

* cell-cell communication;
* molecular signaling;
* immune activation;
* immune suppression;
* functional interaction;
* disease mechanism;
* clinical relevance;
* causality.

For example:

```text
high degree
```

means that a cell has many graph neighbors under the selected spatial rule.

It does not necessarily mean that the cell is functionally influential.

Similarly:

```text
high PageRank
```

means that a cell occupies a prominent position in the defined network.

It does not establish biological control.

MOSAIC therefore distinguishes between:

```text
structural observation
candidate biological interpretation
validated biological conclusion
```

---

# Relationship to Hypothesis-Driven Analysis

Layer 2 initially measures the graph constructed from general spatial relationships such as:

```text
Cell — NEAR_TO — Cell
```

These measurements provide a description of tissue topology.

Later MOSAIC layers will create biologically qualified relationships such as:

```text
TumorCell — TUMOR_CD8_CONTACT — CD8Cell
```

or:

```text
Macrophage — NEAR_TREG — Treg
```

The same graph metrics may then be applied to a hypothesis-specific projection.

This distinction is essential.

A metric calculated on a general proximity graph answers:

> How is the tissue spatially organized?

A metric calculated on a hypothesis-specific graph answers:

> What structure emerges when only relationships relevant to this biological hypothesis are considered?

The graph presented to an algorithm defines the version of biological reality that the algorithm is permitted to analyze.

---

# Recommended Functions

The initial Layer 2 implementation should include:

```python
calculate_graph_qc()
compute_degree()
compute_weighted_degree()
compute_pagerank()
compute_betweenness()
compute_clustering_coefficient()
compute_connected_components()
compute_network_features()
aggregate_metrics_by_cell_type()
aggregate_metrics_by_sample()
```

---

# Proposed Code Structure

```text
mosaic/
└── layer2/
    ├── __init__.py
    ├── qc.py
    ├── centrality.py
    ├── components.py
    └── aggregation.py
```

### `qc.py`

Responsible for:

```python
calculate_graph_qc()
find_isolated_nodes()
find_high_degree_outliers()
edge_distance_summary()
```

### `centrality.py`

Responsible for:

```python
compute_degree()
compute_weighted_degree()
compute_pagerank()
compute_betweenness()
compute_clustering_coefficient()
compute_network_features()
```

### `components.py`

Responsible for:

```python
compute_connected_components()
connected_component_summary()
```

### `aggregation.py`

Responsible for:

```python
aggregate_metrics_by_cell_type()
aggregate_metrics_by_sample()
```

---

# Minimum Example Workflow

```python
import pandas as pd

from mosaic.layer0.validation import validate_cell_table
from mosaic.layer1.radius_graph import build_radius_graph
from mosaic.layer2.qc import calculate_graph_qc
from mosaic.layer2.centrality import compute_network_features
from mosaic.layer2.aggregation import (
    aggregate_metrics_by_cell_type,
    aggregate_metrics_by_sample,
)

cells = pd.read_csv("examples/synthetic_cells.csv")

validated_cells = validate_cell_table(
    cells,
    cell_id_col="cell_id",
    x_col="x",
    y_col="y",
    phenotype_col="cell_type",
    sample_col="sample_id",
)

edges = build_radius_graph(
    validated_cells,
    radius_um=20,
    cell_id_col="cell_id",
    x_col="x",
    y_col="y",
    sample_col="sample_id",
)

qc_report = calculate_graph_qc(
    cells=validated_cells,
    edges=edges,
    cell_id_col="cell_id",
    sample_col="sample_id",
)

cell_metrics = compute_network_features(
    cells=validated_cells,
    edges=edges,
    cell_id_col="cell_id",
    source_col="source",
    target_col="target",
)

cell_type_summary = aggregate_metrics_by_cell_type(
    cell_metrics,
    cell_type_col="cell_type",
    sample_col="sample_id",
)

sample_summary = aggregate_metrics_by_sample(
    cell_metrics,
    sample_col="sample_id",
)

print(qc_report)
print(cell_metrics.head())
print(cell_type_summary.head())
print(sample_summary.head())
```

---

# Outputs

Layer 2 produces:

1. **Graph-quality report**
2. **Cell-level network-feature table**
3. **Cell-type-level network summary**
4. **Sample-level graph summary**
5. **Connected-component assignments**
6. **Flags for possible graph abnormalities**

These outputs become inputs for later MOSAIC layers.

---

# Connection to Later Layers

Layer 2 provides the measurable graph properties required for:

## Layer 3 — Organizational States

Identification of candidate:

* cellular communities;
* hubs;
* interfaces;
* motifs;
* hypothesis-specific connected structures;
* recurrent neighborhood patterns.

## Layer 4 — Tissue Fingerprints

Construction of sample-level organizational signatures using:

* network-feature distributions;
* component structure;
* interaction patterns;
* organizational-state frequencies;
* topology-derived summaries.

## Layer 5 — Multiomic Integration

Integration of network organization with:

* gene expression;
* protein expression;
* mutations;
* pathways;
* clinical outcomes;
* treatment response;
* external biological knowledge.

---

# Design Principles

## 1. Quality before interpretation

The graph must be evaluated before network metrics are interpreted biologically.

## 2. Preserve cell-level results

Aggregated summaries should never replace the underlying cell-level measurements.

## 3. Preserve sample identity

Graph construction, quality assessment, and metric aggregation must remain traceable to the original sample.

## 4. Document graph assumptions

Every result depends on the graph-construction rule.

Parameters such as radius, nearest-neighbor count, weighting method, and filtering criteria must be recorded.

## 5. Separate topology from function

Network position is a structural observation.

Functional meaning requires independent biological evidence.

## 6. Distinguish candidates from validated states

Graph metrics can identify candidate organizational patterns.

They do not independently establish validated biological states.

## 7. Keep metrics interpretable

The initial MOSAIC release prioritizes transparent metrics whose calculations and limitations can be clearly explained.

## 8. Support reproducibility

Identical input data and parameters should produce identical graph-quality reports and network-feature outputs.

---

# Summary

Layer 2 transforms a constructed spatial graph into a quality-controlled and measurable representation of tissue organization.

The workflow is:

```text
Spatial graph
↓
Graph quality assessment
↓
Cell-level network features
↓
Cell-type and sample-level summaries
↓
Candidate structural interpretations
```

Layer 2 does not ask network metrics to decide what biology means.

It evaluates whether the graph is credible, measures the structures that emerge from it, and provides transparent outputs for hypothesis development and independent validation.

**Presence is not architecture.**

**Architecture must first be constructed carefully, evaluated critically, and measured transparently.**
