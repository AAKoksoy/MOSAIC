"""
Network-feature calculation for MOSAIC Layer 2.

Network features are calculated independently for each biological sample.
All cells are retained, including isolated cells, and repeated cell IDs in
different samples are never mixed.
"""

from __future__ import annotations

from collections.abc import Callable

import networkx as nx
import numpy as np
import pandas as pd


def _require_dataframe(value: pd.DataFrame, name: str) -> None:
    """Raise a clear error when an input is not a pandas DataFrame."""

    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame.")


def _require_columns(
    table: pd.DataFrame,
    required_columns: list[str],
    table_name: str,
) -> None:
    """Confirm that a table contains all required columns."""

    missing_columns = [
        column for column in required_columns if column not in table.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {table_name}: "
            + ", ".join(str(column) for column in missing_columns)
        )


def _validate_cells(
    cells: pd.DataFrame,
    *,
    cell_id_col: str,
    sample_col: str,
) -> pd.DataFrame:
    """Validate and standardize the node table."""

    _require_dataframe(cells, "cells")
    _require_columns(cells, [cell_id_col, sample_col], "cells")

    if cells.empty:
        raise ValueError("cells cannot be empty.")

    if cells[[cell_id_col, sample_col]].isna().any().any():
        raise ValueError(
            "Cell identifiers and sample identifiers cannot be missing."
        )

    working_cells = cells[[sample_col, cell_id_col]].copy()
    working_cells[sample_col] = working_cells[sample_col].astype(str)
    working_cells[cell_id_col] = working_cells[cell_id_col].astype(str)

    if (
        working_cells[cell_id_col].str.strip().eq("").any()
        or working_cells[sample_col].str.strip().eq("").any()
    ):
        raise ValueError(
            "Cell identifiers and sample identifiers cannot be empty."
        )

    if working_cells.duplicated(
        subset=[sample_col, cell_id_col]
    ).any():
        raise ValueError(
            "Cell identifiers must be unique within each sample."
        )

    return working_cells


def _validate_edges(
    edges: pd.DataFrame,
    *,
    sample_col: str,
    source_col: str,
    target_col: str,
    weight_col: str | None,
) -> pd.DataFrame:
    """Validate and standardize the edge table."""

    _require_dataframe(edges, "edges")

    required_columns = [sample_col, source_col, target_col]

    if weight_col is not None:
        required_columns.append(weight_col)

    _require_columns(edges, required_columns, "edges")

    working_edges = edges[required_columns].copy()

    if working_edges[
        [sample_col, source_col, target_col]
    ].isna().any().any():
        raise ValueError(
            "Edge sample, source, and target values cannot be missing."
        )

    working_edges[sample_col] = working_edges[sample_col].astype(str)
    working_edges[source_col] = working_edges[source_col].astype(str)
    working_edges[target_col] = working_edges[target_col].astype(str)

    if weight_col is not None:
        working_edges[weight_col] = pd.to_numeric(
            working_edges[weight_col],
            errors="coerce",
        )

        invalid_weight_mask = (
            working_edges[weight_col].isna()
            | ~np.isfinite(working_edges[weight_col])
            | (working_edges[weight_col] < 0)
        )

        if invalid_weight_mask.any():
            raise ValueError(
                f"{weight_col} must contain finite, non-negative numeric values."
            )

    return working_edges


def _build_sample_graphs(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str,
    sample_col: str,
    source_col: str,
    target_col: str,
    weight_col: str | None = None,
) -> list[tuple[str, nx.Graph, list[str]]]:
    """
    Build one undirected NetworkX graph per biological sample.

    Edges referencing nodes absent from the cell table are excluded. These
    integrity problems are reported separately by Layer 2 graph QC.
    """

    working_cells = _validate_cells(
        cells,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
    )

    working_edges = _validate_edges(
        edges,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )

    sample_graphs: list[tuple[str, nx.Graph, list[str]]] = []

    for sample_value, sample_cells in working_cells.groupby(
        sample_col,
        sort=False,
        dropna=False,
    ):
        sample_nodes = sample_cells[cell_id_col].tolist()
        known_nodes = set(sample_nodes)

        sample_edges = working_edges.loc[
            working_edges[sample_col] == sample_value
        ]

        graph = nx.Graph()
        graph.add_nodes_from(sample_nodes)

        if weight_col is None:
            for source_value, target_value in sample_edges[
                [source_col, target_col]
            ].itertuples(index=False, name=None):
                if (
                    source_value in known_nodes
                    and target_value in known_nodes
                ):
                    graph.add_edge(source_value, target_value)
        else:
            for source_value, target_value, weight_value in sample_edges[
                [source_col, target_col, weight_col]
            ].itertuples(index=False, name=None):
                if (
                    source_value in known_nodes
                    and target_value in known_nodes
                ):
                    weight_value = float(weight_value)

                    if graph.has_edge(source_value, target_value):
                        graph[source_value][target_value][
                            weight_col
                        ] += weight_value
                    else:
                        graph.add_edge(
                            source_value,
                            target_value,
                            **{weight_col: weight_value},
                        )

        sample_graphs.append(
            (sample_value, graph, sample_nodes)
        )

    return sample_graphs


def _metric_table(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    metric_name: str,
    metric_calculator: Callable[[nx.Graph], dict[str, float]],
    cell_id_col: str,
    sample_col: str,
    source_col: str,
    target_col: str,
    weight_col: str | None = None,
) -> pd.DataFrame:
    """Apply a NetworkX metric independently to each sample graph."""

    records: list[dict[str, object]] = []

    sample_graphs = _build_sample_graphs(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )

    for sample_value, graph, sample_nodes in sample_graphs:
        metric_values = metric_calculator(graph)

        for node in sample_nodes:
            records.append(
                {
                    sample_col: sample_value,
                    cell_id_col: node,
                    metric_name: float(
                        metric_values.get(node, 0.0)
                    ),
                }
            )

    return pd.DataFrame(
        records,
        columns=[sample_col, cell_id_col, metric_name],
    )


def compute_degree(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
) -> pd.DataFrame:
    """
    Calculate unnormalized degree for every cell.

    Degree is the number of unique neighboring cells. Isolated cells receive
    a degree of zero.
    """

    output = _metric_table(
        cells,
        edges,
        metric_name="degree",
        metric_calculator=lambda graph: dict(graph.degree()),
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    output["degree"] = output["degree"].astype("int64")

    return output


def compute_weighted_degree(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    weight_col: str = "weight",
) -> pd.DataFrame:
    """
    Calculate weighted degree for every cell.

    Weighted degree is the sum of incident edge weights. If duplicate edges
    occur, their weights are added before the metric is calculated.
    """

    return _metric_table(
        cells,
        edges,
        metric_name="weighted_degree",
        metric_calculator=lambda graph: dict(
            graph.degree(weight=weight_col)
        ),
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )


def compute_pagerank(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    weight_col: str = "weight",
    alpha: float = 0.85,
) -> pd.DataFrame:
    """
    Calculate weighted PageRank for every cell.

    PageRank measures a cell's influence based on both its connections and
    the influence of the cells connected to it.
    """

    if not isinstance(
        alpha,
        (int, float, np.integer, np.floating),
    ):
        raise TypeError("alpha must be numeric.")

    alpha = float(alpha)

    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be greater than zero and less than one.")

    def calculate(graph: nx.Graph) -> dict[str, float]:
        if graph.number_of_nodes() == 0:
            return {}

        return nx.pagerank(
            graph,
            alpha=alpha,
            weight=weight_col,
        )

    return _metric_table(
        cells,
        edges,
        metric_name="pagerank",
        metric_calculator=calculate,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )


def compute_betweenness(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    distance_col: str = "distance_um",
    normalized: bool = True,
    approximate_k: int | None = None,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Calculate betweenness centrality for every cell.

    Edge distance is treated as path cost: shorter spatial edges contribute
    shorter paths. ``approximate_k`` may be used to sample source nodes for
    faster calculation on large graphs.
    """

    if not isinstance(normalized, bool):
        raise TypeError("normalized must be a boolean.")

    if approximate_k is not None:
        if (
            not isinstance(approximate_k, (int, np.integer))
            or isinstance(approximate_k, bool)
        ):
            raise TypeError("approximate_k must be an integer or None.")

        if approximate_k < 1:
            raise ValueError(
                "approximate_k must be greater than or equal to one."
            )

    def calculate(graph: nx.Graph) -> dict[str, float]:
        node_count = graph.number_of_nodes()

        if node_count == 0:
            return {}

        sample_size = approximate_k

        if sample_size is not None:
            sample_size = min(sample_size, node_count)

        return nx.betweenness_centrality(
            graph,
            k=sample_size,
            normalized=normalized,
            weight=distance_col,
            seed=random_seed,
        )

    return _metric_table(
        cells,
        edges,
        metric_name="betweenness",
        metric_calculator=calculate,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=distance_col,
    )


def compute_clustering_coefficient(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
) -> pd.DataFrame:
    """
    Calculate the local clustering coefficient for every cell.

    The coefficient measures how frequently a cell's neighbors are also
    connected to one another. Cells with fewer than two neighbors receive
    a value of zero.
    """

    return _metric_table(
        cells,
        edges,
        metric_name="clustering_coefficient",
        metric_calculator=nx.clustering,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )


def compute_network_features(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    weight_col: str = "weight",
    distance_col: str = "distance_um",
    pagerank_alpha: float = 0.85,
    approximate_betweenness_k: int | None = None,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Calculate the complete MOSAIC Layer 2 cell-level feature table.

    Returned features are:

    - degree,
    - weighted degree,
    - PageRank,
    - betweenness centrality,
    - local clustering coefficient.
    """

    degree = compute_degree(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    weighted_degree = compute_weighted_degree(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
    )

    pagerank = compute_pagerank(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        weight_col=weight_col,
        alpha=pagerank_alpha,
    )

    betweenness = compute_betweenness(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        distance_col=distance_col,
        approximate_k=approximate_betweenness_k,
        random_seed=random_seed,
    )

    clustering = compute_clustering_coefficient(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    key_columns = [sample_col, cell_id_col]

    feature_table = degree.merge(
        weighted_degree,
        on=key_columns,
        how="left",
        validate="one_to_one",
    )

    for metric_table in [
        pagerank,
        betweenness,
        clustering,
    ]:
        feature_table = feature_table.merge(
            metric_table,
            on=key_columns,
            how="left",
            validate="one_to_one",
        )

    metric_columns = [
        "degree",
        "weighted_degree",
        "pagerank",
        "betweenness",
        "clustering_coefficient",
    ]

    feature_table[metric_columns] = feature_table[
        metric_columns
    ].fillna(0.0)

    feature_table["degree"] = feature_table["degree"].astype(
        "int64"
    )

    return feature_table[
        key_columns + metric_columns
    ]
