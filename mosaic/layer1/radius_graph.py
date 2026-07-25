"""
Radius-based spatial graph construction for MOSAIC Layer 1.

This module converts validated cell-level observations into an
undirected spatial edge table. Cells are connected when the Euclidean
distance between them is less than or equal to a user-defined radius.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from mosaic.layer0.validation import validate_cell_table


REQUIRED_EDGE_COLUMNS = (
    "source",
    "target",
    "distance_um",
    "sample_id",
    "relationship_type",
)


def build_radius_graph(
    cells: pd.DataFrame,
    *,
    radius_um: float = 20.0,
    cell_id_col: str = "cell_id",
    x_col: str = "x",
    y_col: str = "y",
    phenotype_col: str = "cell_type",
    sample_col: str = "sample_id",
    relationship_type: str = "NEAR_TO",
    include_edge_weight: bool = False,
    weight_epsilon: float = 1e-9,
) -> pd.DataFrame:
    """
    Construct an undirected spatial radius graph.

    Cells are connected when they belong to the same sample and their
    Euclidean distance is less than or equal to ``radius_um``.

    Parameters
    ----------
    cells:
        Cell-level table containing one row per biological entity.
    radius_um:
        Maximum Euclidean distance permitted between connected cells.
        The value must be greater than zero.
    cell_id_col:
        Column containing the cell identifier.
    x_col:
        Column containing the x coordinate.
    y_col:
        Column containing the y coordinate.
    phenotype_col:
        Column containing the cell phenotype or cell type.
    sample_col:
        Column containing the sample identifier.
    relationship_type:
        Label assigned to every spatial relationship.
    include_edge_weight:
        Whether to include an inverse-distance edge weight.
    weight_epsilon:
        Small positive value added to distance when calculating
        inverse-distance weights.

    Returns
    -------
    pandas.DataFrame
        Spatial edge table containing one row per undirected cell pair.

    Notes
    -----
    Graphs are constructed independently within each sample. Cells from
    different samples are never connected.

    Each undirected cell pair is stored once. Self-connections are excluded.

    The default edge weight, when requested, is:

    ``1 / (distance_um + weight_epsilon)``

    This weighting function gives stronger weights to closer cells.
    """

    if not isinstance(radius_um, (int, float, np.integer, np.floating)):
        raise TypeError("radius_um must be numeric.")

    radius_um = float(radius_um)

    if not np.isfinite(radius_um) or radius_um <= 0:
        raise ValueError("radius_um must be a finite value greater than zero.")

    if not isinstance(relationship_type, str):
        raise TypeError("relationship_type must be a string.")

    relationship_type = relationship_type.strip()

    if relationship_type == "":
        raise ValueError("relationship_type cannot be empty.")

    if not isinstance(weight_epsilon, (int, float, np.integer, np.floating)):
        raise TypeError("weight_epsilon must be numeric.")

    weight_epsilon = float(weight_epsilon)

    if not np.isfinite(weight_epsilon) or weight_epsilon <= 0:
        raise ValueError(
            "weight_epsilon must be a finite value greater than zero."
        )

    validated_cells = validate_cell_table(
        cells,
        cell_id_col=cell_id_col,
        x_col=x_col,
        y_col=y_col,
        phenotype_col=phenotype_col,
        sample_col=sample_col,
    )

    edge_records: list[dict[str, Any]] = []

    grouped_cells = validated_cells.groupby(
        sample_col,
        sort=False,
        dropna=False,
    )

    for sample_id, sample_cells in grouped_cells:
        if len(sample_cells) < 2:
            continue

        sample_cells = sample_cells.reset_index(drop=True)

        coordinates = sample_cells[[x_col, y_col]].to_numpy(
            dtype=float,
            copy=True,
        )

        cell_ids = sample_cells[cell_id_col].to_numpy(copy=True)

        spatial_tree = cKDTree(coordinates)

        neighbor_pairs = spatial_tree.query_pairs(
            r=radius_um,
            output_type="set",
        )

        for first_index, second_index in sorted(neighbor_pairs):
            source_id = cell_ids[first_index]
            target_id = cell_ids[second_index]

            source_coordinates = coordinates[first_index]
            target_coordinates = coordinates[second_index]

            distance_um = float(
                np.linalg.norm(source_coordinates - target_coordinates)
            )

            edge_record: dict[str, Any] = {
                "source": source_id,
                "target": target_id,
                "distance_um": distance_um,
                "sample_id": sample_id,
                "relationship_type": relationship_type,
            }

            if include_edge_weight:
                edge_record["edge_weight"] = float(
                    1.0 / (distance_um + weight_epsilon)
                )

            edge_records.append(edge_record)

    output_columns = list(REQUIRED_EDGE_COLUMNS)

    if include_edge_weight:
        output_columns.append("edge_weight")

    edges = pd.DataFrame(
        edge_records,
        columns=output_columns,
    )

    return edges


def validate_edge_table(
    edges: pd.DataFrame,
    *,
    source_col: str = "source",
    target_col: str = "target",
    distance_col: str = "distance_um",
    sample_col: str = "sample_id",
    relationship_col: str = "relationship_type",
    radius_um: float | None = None,
    numerical_tolerance: float = 1e-9,
) -> pd.DataFrame:
    """
    Validate a MOSAIC spatial edge table.

    Parameters
    ----------
    edges:
        Spatial edge table.
    source_col:
        Column containing the source cell identifier.
    target_col:
        Column containing the target cell identifier.
    distance_col:
        Column containing Euclidean edge distance.
    sample_col:
        Column containing the sample identifier.
    relationship_col:
        Column containing the relationship label.
    radius_um:
        Optional maximum permitted edge distance.
    numerical_tolerance:
        Tolerance used when comparing distances with ``radius_um``.

    Returns
    -------
    pandas.DataFrame
        A validated copy of the edge table.

    Raises
    ------
    TypeError
        If the input is not a pandas DataFrame.
    ValueError
        If required columns are missing or invalid edges are found.
    """

    if not isinstance(edges, pd.DataFrame):
        raise TypeError("edges must be a pandas DataFrame.")

    required_columns = [
        source_col,
        target_col,
        distance_col,
        sample_col,
        relationship_col,
    ]

    missing_columns = [
        column for column in required_columns if column not in edges.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required edge columns: "
            + ", ".join(str(column) for column in missing_columns)
        )

    validated = edges.copy()

    if validated.empty:
        return validated

    missing_required_values = validated[required_columns].isna().sum()

    columns_with_missing_values = {
        column: int(count)
        for column, count in missing_required_values.items()
        if count > 0
    }

    if columns_with_missing_values:
        details = ", ".join(
            f"{column}={count}"
            for column, count in columns_with_missing_values.items()
        )

        raise ValueError(
            "The edge table contains missing required values: " + details
        )

    numeric_distances = pd.to_numeric(
        validated[distance_col],
        errors="coerce",
    )

    invalid_distance_count = int(numeric_distances.isna().sum())

    if invalid_distance_count > 0:
        raise ValueError(
            f"Column '{distance_col}' contains "
            f"{invalid_distance_count} nonnumeric value(s)."
        )

    if not np.isfinite(numeric_distances.to_numpy(dtype=float)).all():
        raise ValueError(
            f"Column '{distance_col}' contains non-finite values."
        )

    if (numeric_distances < 0).any():
        raise ValueError("Edge distances cannot be negative.")

    validated[distance_col] = numeric_distances.astype(float)

    self_loop_mask = (
        validated[source_col].astype(str)
        == validated[target_col].astype(str)
    )

    if self_loop_mask.any():
        self_loop_count = int(self_loop_mask.sum())

        raise ValueError(
            f"Found {self_loop_count} self-loop edge(s)."
        )

    duplicate_mask = _find_duplicate_undirected_edges(
        validated,
        source_col=source_col,
        target_col=target_col,
        sample_col=sample_col,
        relationship_col=relationship_col,
    )

    if duplicate_mask.any():
        duplicate_count = int(duplicate_mask.sum())

        raise ValueError(
            f"Found {duplicate_count} duplicated undirected edge row(s)."
        )

    if radius_um is not None:
        if not isinstance(
            radius_um,
            (int, float, np.integer, np.floating),
        ):
            raise TypeError("radius_um must be numeric when provided.")

        radius_um = float(radius_um)

        if not np.isfinite(radius_um) or radius_um <= 0:
            raise ValueError(
                "radius_um must be a finite value greater than zero."
            )

        excessive_distance_mask = (
            validated[distance_col]
            > radius_um + numerical_tolerance
        )

        if excessive_distance_mask.any():
            excessive_count = int(excessive_distance_mask.sum())

            raise ValueError(
                f"Found {excessive_count} edge(s) with distance greater "
                f"than the specified radius of {radius_um}."
            )

    return validated


def summarize_spatial_graph(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
    distance_col: str = "distance_um",
) -> dict[str, Any]:
    """
    Generate a basic summary of a spatial graph.

    Parameters
    ----------
    cells:
        Validated cell table.
    edges:
        Spatial edge table.
    cell_id_col:
        Column containing the cell identifier.
    sample_col:
        Column containing the sample identifier.
    source_col:
        Column containing source cell identifiers.
    target_col:
        Column containing target cell identifiers.
    distance_col:
        Column containing edge distances.

    Returns
    -------
    dict
        Summary statistics describing the spatial graph.
    """

    if not isinstance(cells, pd.DataFrame):
        raise TypeError("cells must be a pandas DataFrame.")

    if not isinstance(edges, pd.DataFrame):
        raise TypeError("edges must be a pandas DataFrame.")

    required_cell_columns = [cell_id_col, sample_col]

    missing_cell_columns = [
        column
        for column in required_cell_columns
        if column not in cells.columns
    ]

    if missing_cell_columns:
        raise ValueError(
            "Missing required cell columns: "
            + ", ".join(str(column) for column in missing_cell_columns)
        )

    validated_edges = validate_edge_table(
        edges,
        source_col=source_col,
        target_col=target_col,
        distance_col=distance_col,
        sample_col=sample_col,
    )

    number_of_nodes = int(len(cells))
    number_of_edges = int(len(validated_edges))

    if number_of_nodes > 1:
        graph_density = float(
            (2.0 * number_of_edges)
            / (number_of_nodes * (number_of_nodes - 1))
        )
    else:
        graph_density = 0.0

    degree_counts = _calculate_degree_counts(
        cells=cells,
        edges=validated_edges,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
        source_col=source_col,
        target_col=target_col,
    )

    isolated_nodes = int((degree_counts == 0).sum())

    if number_of_nodes > 0:
        fraction_isolated = float(isolated_nodes / number_of_nodes)
    else:
        fraction_isolated = 0.0

    distance_values = validated_edges[distance_col]

    summary = {
        "number_of_nodes": number_of_nodes,
        "number_of_edges": number_of_edges,
        "number_of_samples": int(
            cells[sample_col].nunique(dropna=True)
        ),
        "graph_density": graph_density,
        "number_of_isolated_nodes": isolated_nodes,
        "fraction_isolated": fraction_isolated,
        "minimum_degree": _safe_numeric_min(degree_counts),
        "maximum_degree": _safe_numeric_max(degree_counts),
        "mean_degree": _safe_numeric_mean(degree_counts),
        "median_degree": _safe_numeric_median(degree_counts),
        "minimum_edge_distance": _safe_numeric_min(distance_values),
        "maximum_edge_distance": _safe_numeric_max(distance_values),
        "mean_edge_distance": _safe_numeric_mean(distance_values),
        "median_edge_distance": _safe_numeric_median(distance_values),
    }

    return summary


def _find_duplicate_undirected_edges(
    edges: pd.DataFrame,
    *,
    source_col: str,
    target_col: str,
    sample_col: str,
    relationship_col: str,
) -> pd.Series:
    """
    Identify duplicate undirected edges.

    Source-target order is ignored when determining duplication.
    """

    source_values = edges[source_col].astype(str)
    target_values = edges[target_col].astype(str)

    canonical_source = np.minimum(source_values, target_values)
    canonical_target = np.maximum(source_values, target_values)

    duplicate_frame = pd.DataFrame(
        {
            "sample": edges[sample_col].astype(str),
            "relationship": edges[relationship_col].astype(str),
            "canonical_source": canonical_source,
            "canonical_target": canonical_target,
        },
        index=edges.index,
    )

    return duplicate_frame.duplicated(
        subset=[
            "sample",
            "relationship",
            "canonical_source",
            "canonical_target",
        ],
        keep=False,
    )


def _calculate_degree_counts(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str,
    sample_col: str,
    source_col: str,
    target_col: str,
) -> pd.Series:
    """
    Calculate degree while preserving sample-specific cell identity.
    """

    node_index = pd.MultiIndex.from_frame(
        cells[[sample_col, cell_id_col]].astype(str)
    )

    degree_counts = pd.Series(
        0,
        index=node_index,
        dtype=int,
    )

    if edges.empty:
        return degree_counts

    source_index = pd.MultiIndex.from_arrays(
        [
            edges[sample_col].astype(str),
            edges[source_col].astype(str),
        ]
    )

    target_index = pd.MultiIndex.from_arrays(
        [
            edges[sample_col].astype(str),
            edges[target_col].astype(str),
        ]
    )

    source_counts = pd.Series(
        1,
        index=source_index,
    ).groupby(level=[0, 1]).sum()

    target_counts = pd.Series(
        1,
        index=target_index,
    ).groupby(level=[0, 1]).sum()

    total_counts = source_counts.add(
        target_counts,
        fill_value=0,
    ).astype(int)

    degree_counts.loc[
        degree_counts.index.intersection(total_counts.index)
    ] = total_counts.loc[
        degree_counts.index.intersection(total_counts.index)
    ]

    return degree_counts


def _safe_numeric_min(values: pd.Series) -> float | int | None:
    """Return the minimum numeric value or None when unavailable."""

    numeric_values = pd.to_numeric(values, errors="coerce").dropna()

    if numeric_values.empty:
        return None

    return float(numeric_values.min())


def _safe_numeric_max(values: pd.Series) -> float | int | None:
    """Return the maximum numeric value or None when unavailable."""

    numeric_values = pd.to_numeric(values, errors="coerce").dropna()

    if numeric_values.empty:
        return None

    return float(numeric_values.max())


def _safe_numeric_mean(values: pd.Series) -> float | None:
    """Return the numeric mean or None when unavailable."""

    numeric_values = pd.to_numeric(values, errors="coerce").dropna()

    if numeric_values.empty:
        return None

    return float(numeric_values.mean())


def _safe_numeric_median(values: pd.Series) -> float | None:
    """Return the numeric median or None when unavailable."""

    numeric_values = pd.to_numeric(values, errors="coerce").dropna()

    if numeric_values.empty:
        return None

    return float(numeric_values.median())
