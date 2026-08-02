"""
Graph quality-control utilities for MOSAIC Layer 2.

This module evaluates spatial graph integrity, identifies isolated cells
and degree outliers, and summarizes graph quality independently for each
biological sample.
"""

from __future__ import annotations

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


def _calculate_degree_table(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str,
    sample_col: str,
    source_col: str,
    target_col: str,
) -> pd.DataFrame:
    """
    Calculate undirected degree while preserving sample boundaries.

    The MultiIndex level names are explicitly standardized so source and
    target counts align correctly on ``(sample_id, node)``.
    """

    node_index = pd.MultiIndex.from_arrays(
        [
            cells[sample_col].astype(str),
            cells[cell_id_col].astype(str),
        ],
        names=[sample_col, "node"],
    )

    degree_counts = pd.Series(
        0,
        index=node_index,
        dtype="int64",
        name="degree",
    )

    if not edges.empty:
        source_index = pd.MultiIndex.from_arrays(
            [
                edges[sample_col].astype(str),
                edges[source_col].astype(str),
            ],
            names=[sample_col, "node"],
        )

        target_index = pd.MultiIndex.from_arrays(
            [
                edges[sample_col].astype(str),
                edges[target_col].astype(str),
            ],
            names=[sample_col, "node"],
        )

        source_counts = pd.Series(
            1,
            index=source_index,
            dtype="int64",
        ).groupby(level=[0, 1]).sum()

        target_counts = pd.Series(
            1,
            index=target_index,
            dtype="int64",
        ).groupby(level=[0, 1]).sum()

        combined_counts = source_counts.add(
            target_counts,
            fill_value=0,
        ).astype("int64")

        matching_nodes = degree_counts.index.intersection(
            combined_counts.index
        )

        degree_counts.loc[matching_nodes] = combined_counts.loc[
            matching_nodes
        ]

    degree_table = degree_counts.rename("degree").reset_index()

    degree_table[sample_col] = degree_table[sample_col].astype(str)
    degree_table["node"] = degree_table["node"].astype(str)

    return degree_table


def find_isolated_nodes(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
) -> pd.DataFrame:
    """
    Return cells with degree zero.

    Isolation is evaluated within each sample, so identical cell IDs may
    safely appear in different samples.

    Parameters
    ----------
    cells:
        Cell-level table containing all graph nodes.
    edges:
        Spatial edge table.
    cell_id_col:
        Column containing cell identifiers.
    sample_col:
        Column containing sample identifiers.
    source_col:
        Edge source column.
    target_col:
        Edge target column.

    Returns
    -------
    pandas.DataFrame
        One row per isolated cell, including ``sample_id``, ``cell_id``,
        and ``degree``.
    """

    _require_dataframe(cells, "cells")
    _require_dataframe(edges, "edges")

    _require_columns(
        cells,
        [cell_id_col, sample_col],
        "cells",
    )

    _require_columns(
        edges,
        [source_col, target_col, sample_col],
        "edges",
    )

    degree_table = _calculate_degree_table(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    isolated_nodes = degree_table.loc[
        degree_table["degree"] == 0
    ].copy()

    isolated_nodes = isolated_nodes.rename(
        columns={"node": cell_id_col}
    )

    return isolated_nodes.reset_index(drop=True)


def find_high_degree_outliers(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    iqr_multiplier: float = 1.5,
) -> pd.DataFrame:
    """
    Identify unusually high-degree cells using a sample-wise IQR rule.

    A cell is classified as an outlier when:

    ``degree > Q3 + iqr_multiplier * IQR``

    Parameters
    ----------
    cells:
        Cell-level table containing all graph nodes.
    edges:
        Spatial edge table.
    cell_id_col:
        Column containing cell identifiers.
    sample_col:
        Column containing sample identifiers.
    source_col:
        Edge source column.
    target_col:
        Edge target column.
    iqr_multiplier:
        Multiplier applied to the interquartile range. The conventional
        value is 1.5.

    Returns
    -------
    pandas.DataFrame
        High-degree cells with their sample-specific thresholds.
    """

    _require_dataframe(cells, "cells")
    _require_dataframe(edges, "edges")

    _require_columns(
        cells,
        [cell_id_col, sample_col],
        "cells",
    )

    _require_columns(
        edges,
        [source_col, target_col, sample_col],
        "edges",
    )

    if not isinstance(
        iqr_multiplier,
        (int, float, np.integer, np.floating),
    ):
        raise TypeError("iqr_multiplier must be numeric.")

    iqr_multiplier = float(iqr_multiplier)

    if not np.isfinite(iqr_multiplier) or iqr_multiplier < 0:
        raise ValueError(
            "iqr_multiplier must be finite and greater than or equal to zero."
        )

    degree_table = _calculate_degree_table(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    output_tables: list[pd.DataFrame] = []

    for _, sample_degrees in degree_table.groupby(
        sample_col,
        sort=False,
        dropna=False,
    ):
        sample_degrees = sample_degrees.copy()

        first_quartile = float(
            sample_degrees["degree"].quantile(0.25)
        )
        third_quartile = float(
            sample_degrees["degree"].quantile(0.75)
        )

        interquartile_range = third_quartile - first_quartile
        threshold = third_quartile + (
            iqr_multiplier * interquartile_range
        )

        sample_degrees["degree_outlier_threshold"] = threshold

        sample_outliers = sample_degrees.loc[
            sample_degrees["degree"] > threshold
        ].copy()

        output_tables.append(sample_outliers)

    output_columns = [
        sample_col,
        cell_id_col,
        "degree",
        "degree_outlier_threshold",
    ]

    if not output_tables:
        return pd.DataFrame(columns=output_columns)

    outliers = pd.concat(
        output_tables,
        ignore_index=True,
    )

    outliers = outliers.rename(columns={"node": cell_id_col})

    return outliers[output_columns]


def calculate_graph_qc(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    distance_col: str = "distance_um",
    iqr_multiplier: float = 1.5,
) -> pd.DataFrame:
    """
    Calculate sample-level spatial graph quality metrics.

    The report includes:

    - node and edge counts,
    - graph density,
    - isolated cells,
    - high-degree outliers,
    - self-loops,
    - duplicate undirected edges,
    - edges whose nodes are missing from the cell table,
    - cross-sample violations,
    - minimum, mean, median, and maximum edge distance.

    Parameters
    ----------
    cells:
        Cell-level table containing all graph nodes.
    edges:
        Spatial edge table.
    cell_id_col:
        Column containing cell identifiers.
    sample_col:
        Column containing sample identifiers.
    source_col:
        Edge source column.
    target_col:
        Edge target column.
    distance_col:
        Column containing spatial edge distances.
    iqr_multiplier:
        Multiplier used for sample-wise high-degree outlier detection.

    Returns
    -------
    pandas.DataFrame
        One graph quality-control record per sample.

    Notes
    -----
    Graph density is calculated for a simple undirected graph:

    ``2E / (N * (N - 1))``

    Density is reported as zero for samples containing fewer than two
    nodes.
    """

    _require_dataframe(cells, "cells")
    _require_dataframe(edges, "edges")

    _require_columns(
        cells,
        [cell_id_col, sample_col],
        "cells",
    )

    _require_columns(
        edges,
        [
            source_col,
            target_col,
            sample_col,
            distance_col,
        ],
        "edges",
    )

    if cells.empty:
        raise ValueError("cells cannot be empty.")

    if cells[[cell_id_col, sample_col]].isna().any().any():
        raise ValueError(
            "Cell identifiers and sample identifiers cannot be missing."
        )

    working_cells = cells.copy()
    working_edges = edges.copy()

    working_cells["_qc_sample"] = (
        working_cells[sample_col].astype(str)
    )
    working_cells["_qc_cell"] = (
        working_cells[cell_id_col].astype(str)
    )

    working_edges["_qc_sample"] = (
        working_edges[sample_col].astype(str)
    )
    working_edges["_qc_source"] = (
        working_edges[source_col].astype(str)
    )
    working_edges["_qc_target"] = (
        working_edges[target_col].astype(str)
    )

    numeric_distances = pd.to_numeric(
        working_edges[distance_col],
        errors="coerce",
    )

    invalid_distance_mask = (
        numeric_distances.isna()
        | ~np.isfinite(numeric_distances)
        | (numeric_distances < 0)
    )

    working_edges["_qc_distance"] = numeric_distances

    degree_table = _calculate_degree_table(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    isolated_nodes = find_isolated_nodes(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    high_degree_outliers = find_high_degree_outliers(
        cells,
        edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
        iqr_multiplier=iqr_multiplier,
    )

    known_nodes = set(
        zip(
            working_cells["_qc_sample"],
            working_cells["_qc_cell"],
        )
    )

    source_keys = list(
        zip(
            working_edges["_qc_sample"],
            working_edges["_qc_source"],
        )
    )

    target_keys = list(
        zip(
            working_edges["_qc_sample"],
            working_edges["_qc_target"],
        )
    )

    missing_source_mask = pd.Series(
        [key not in known_nodes for key in source_keys],
        index=working_edges.index,
        dtype=bool,
    )

    missing_target_mask = pd.Series(
        [key not in known_nodes for key in target_keys],
        index=working_edges.index,
        dtype=bool,
    )

    working_edges["_qc_missing_node"] = (
        missing_source_mask | missing_target_mask
    )

    cell_samples: dict[str, set[str]] = {}

    for sample_value, cell_value in zip(
        working_cells["_qc_sample"],
        working_cells["_qc_cell"],
    ):
        cell_samples.setdefault(cell_value, set()).add(sample_value)

    cross_sample_flags: list[bool] = []

    for edge_sample, source_value, target_value in zip(
        working_edges["_qc_sample"],
        working_edges["_qc_source"],
        working_edges["_qc_target"],
    ):
        source_samples = cell_samples.get(source_value, set())
        target_samples = cell_samples.get(target_value, set())

        cross_sample_flags.append(
            (
                bool(source_samples)
                and edge_sample not in source_samples
            )
            or (
                bool(target_samples)
                and edge_sample not in target_samples
            )
        )

    working_edges["_qc_cross_sample"] = pd.Series(
        cross_sample_flags,
        index=working_edges.index,
        dtype=bool,
    )

    working_edges["_qc_self_loop"] = (
        working_edges["_qc_source"]
        == working_edges["_qc_target"]
    )

    working_edges["_qc_node_low"] = working_edges[
        ["_qc_source", "_qc_target"]
    ].min(axis=1)

    working_edges["_qc_node_high"] = working_edges[
        ["_qc_source", "_qc_target"]
    ].max(axis=1)

    working_edges["_qc_duplicate"] = working_edges.duplicated(
        subset=[
            "_qc_sample",
            "_qc_node_low",
            "_qc_node_high",
        ],
        keep="first",
    )

    working_edges["_qc_invalid_distance"] = invalid_distance_mask

    qc_records: list[dict[str, object]] = []

    for sample_value, sample_cells in working_cells.groupby(
        "_qc_sample",
        sort=False,
        dropna=False,
    ):
        sample_edges = working_edges.loc[
            working_edges["_qc_sample"] == sample_value
        ].copy()

        sample_degree_table = degree_table.loc[
            degree_table[sample_col] == sample_value
        ]

        sample_isolated = isolated_nodes.loc[
            isolated_nodes[sample_col] == sample_value
        ]

        sample_outliers = high_degree_outliers.loc[
            high_degree_outliers[sample_col] == sample_value
        ]

        node_count = int(len(sample_cells))
        edge_count = int(len(sample_edges))

        if node_count < 2:
            graph_density = 0.0
        else:
            graph_density = float(
                (2.0 * edge_count)
                / (node_count * (node_count - 1))
            )

        valid_distances = sample_edges.loc[
            ~sample_edges["_qc_invalid_distance"],
            "_qc_distance",
        ]

        if valid_distances.empty:
            minimum_distance = np.nan
            mean_distance = np.nan
            median_distance = np.nan
            maximum_distance = np.nan
        else:
            minimum_distance = float(valid_distances.min())
            mean_distance = float(valid_distances.mean())
            median_distance = float(valid_distances.median())
            maximum_distance = float(valid_distances.max())

        qc_records.append(
            {
                sample_col: sample_value,
                "node_count": node_count,
                "edge_count": edge_count,
                "graph_density": graph_density,
                "mean_degree": float(
                    sample_degree_table["degree"].mean()
                ),
                "median_degree": float(
                    sample_degree_table["degree"].median()
                ),
                "maximum_degree": int(
                    sample_degree_table["degree"].max()
                ),
                "isolated_node_count": int(len(sample_isolated)),
                "isolated_node_fraction": float(
                    len(sample_isolated) / node_count
                ),
                "high_degree_outlier_count": int(
                    len(sample_outliers)
                ),
                "self_loop_count": int(
                    sample_edges["_qc_self_loop"].sum()
                ),
                "duplicate_edge_count": int(
                    sample_edges["_qc_duplicate"].sum()
                ),
                "missing_node_edge_count": int(
                    sample_edges["_qc_missing_node"].sum()
                ),
                "cross_sample_violation_count": int(
                    sample_edges["_qc_cross_sample"].sum()
                ),
                "invalid_distance_count": int(
                    sample_edges["_qc_invalid_distance"].sum()
                ),
                "minimum_distance_um": minimum_distance,
                "mean_distance_um": mean_distance,
                "median_distance_um": median_distance,
                "maximum_distance_um": maximum_distance,
            }
        )

    return pd.DataFrame(qc_records)
