"""
Feature aggregation utilities for MOSAIC Layer 2.

This module converts cell-level network features into interpretable summaries
for cell types and biological samples.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


DEFAULT_METRIC_COLUMNS = [
    "degree",
    "weighted_degree",
    "pagerank",
    "betweenness",
    "clustering_coefficient",
]


def _require_dataframe(
    value: pd.DataFrame,
    name: str,
) -> None:
    """Raise a clear error when an input is not a pandas DataFrame."""

    if not isinstance(value, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame.")


def _require_columns(
    table: pd.DataFrame,
    required_columns: Sequence[str],
    table_name: str,
) -> None:
    """Confirm that a table contains all required columns."""

    missing_columns = [
        column
        for column in required_columns
        if column not in table.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns in {table_name}: "
            + ", ".join(str(column) for column in missing_columns)
        )


def _validate_metric_columns(
    features: pd.DataFrame,
    metric_columns: Sequence[str] | None,
) -> list[str]:
    """Validate requested network metrics and return them as a list."""

    if metric_columns is None:
        selected_metrics = [
            column
            for column in DEFAULT_METRIC_COLUMNS
            if column in features.columns
        ]
    else:
        if isinstance(metric_columns, str):
            raise TypeError(
                "metric_columns must be a sequence of column names, "
                "not a single string."
            )

        selected_metrics = list(metric_columns)

    if not selected_metrics:
        raise ValueError(
            "No network metric columns were selected or found."
        )

    if len(selected_metrics) != len(set(selected_metrics)):
        raise ValueError(
            "metric_columns cannot contain duplicate column names."
        )

    _require_columns(
        features,
        selected_metrics,
        "features",
    )

    return selected_metrics


def _prepare_features(
    features: pd.DataFrame,
    *,
    metric_columns: Sequence[str] | None,
    cell_id_col: str,
    sample_col: str,
) -> tuple[pd.DataFrame, list[str]]:
    """Validate and standardize a cell-level network-feature table."""

    _require_dataframe(features, "features")
    _require_columns(
        features,
        [sample_col, cell_id_col],
        "features",
    )

    if features.empty:
        raise ValueError("features cannot be empty.")

    if features[[sample_col, cell_id_col]].isna().any().any():
        raise ValueError(
            "Feature sample identifiers and cell identifiers "
            "cannot be missing."
        )

    selected_metrics = _validate_metric_columns(
        features,
        metric_columns,
    )

    working_features = features[
        [sample_col, cell_id_col] + selected_metrics
    ].copy()

    working_features[sample_col] = (
        working_features[sample_col].astype(str)
    )
    working_features[cell_id_col] = (
        working_features[cell_id_col].astype(str)
    )

    if working_features.duplicated(
        subset=[sample_col, cell_id_col]
    ).any():
        raise ValueError(
            "The feature table must contain one row per cell "
            "within each sample."
        )

    for metric in selected_metrics:
        working_features[metric] = pd.to_numeric(
            working_features[metric],
            errors="coerce",
        )

        invalid_metric_mask = (
            working_features[metric].isna()
            | ~np.isfinite(working_features[metric])
        )

        if invalid_metric_mask.any():
            raise ValueError(
                f"{metric} must contain finite numeric values."
            )

    return working_features, selected_metrics


def _flatten_aggregated_columns(
    table: pd.DataFrame,
) -> pd.DataFrame:
    """Convert pandas aggregation MultiIndex columns to flat names."""

    flattened_columns: list[str] = []

    for column in table.columns:
        if isinstance(column, tuple):
            parts = [
                str(part)
                for part in column
                if str(part) not in {"", "None"}
            ]
            flattened_columns.append("_".join(parts))
        else:
            flattened_columns.append(str(column))

    output = table.copy()
    output.columns = flattened_columns

    return output


def _aggregate_metrics(
    table: pd.DataFrame,
    *,
    group_columns: list[str],
    metric_columns: list[str],
) -> pd.DataFrame:
    """Calculate standardized descriptive statistics."""

    aggregated = (
        table.groupby(
            group_columns,
            sort=False,
            dropna=False,
        )[metric_columns]
        .agg(["mean", "median", "std", "min", "max"])
        .reset_index()
    )

    aggregated = _flatten_aggregated_columns(aggregated)

    standard_deviation_columns = [
        f"{metric}_std"
        for metric in metric_columns
    ]

    aggregated[standard_deviation_columns] = aggregated[
        standard_deviation_columns
    ].fillna(0.0)

    return aggregated


def aggregate_metrics_by_cell_type(
    features: pd.DataFrame,
    cells: pd.DataFrame,
    *,
    metric_columns: Sequence[str] | None = None,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    cell_type_col: str = "cell_type",
) -> pd.DataFrame:
    """
    Aggregate cell-level network metrics by sample and cell type.

    The output describes the network position of each cell population within
    each biological sample. Cell counts and within-sample cell fractions are
    included alongside mean, median, standard deviation, minimum, and maximum
    values for every selected metric.

    Parameters
    ----------
    features:
        Cell-level output from ``compute_network_features``.
    cells:
        Layer 0 cell table containing cell-type annotations.
    metric_columns:
        Network metrics to aggregate. If omitted, all available default
        Layer 2 metrics are used.
    cell_id_col:
        Column containing cell identifiers.
    sample_col:
        Column containing sample identifiers.
    cell_type_col:
        Column containing cell-type or phenotype annotations.

    Returns
    -------
    pandas.DataFrame
        One row per sample and cell type.
    """

    working_features, selected_metrics = _prepare_features(
        features,
        metric_columns=metric_columns,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
    )

    _require_dataframe(cells, "cells")
    _require_columns(
        cells,
        [sample_col, cell_id_col, cell_type_col],
        "cells",
    )

    if cells.empty:
        raise ValueError("cells cannot be empty.")

    working_cells = cells[
        [sample_col, cell_id_col, cell_type_col]
    ].copy()

    if working_cells[
        [sample_col, cell_id_col, cell_type_col]
    ].isna().any().any():
        raise ValueError(
            "Cell sample identifiers, cell identifiers, and "
            "cell-type annotations cannot be missing."
        )

    working_cells[sample_col] = (
        working_cells[sample_col].astype(str)
    )
    working_cells[cell_id_col] = (
        working_cells[cell_id_col].astype(str)
    )
    working_cells[cell_type_col] = (
        working_cells[cell_type_col].astype(str)
    )

    if working_cells.duplicated(
        subset=[sample_col, cell_id_col]
    ).any():
        raise ValueError(
            "Cell identifiers must be unique within each sample."
        )

    merged = working_features.merge(
        working_cells,
        on=[sample_col, cell_id_col],
        how="left",
        validate="one_to_one",
        indicator=True,
    )

    unmatched_feature_mask = merged["_merge"] != "both"

    if unmatched_feature_mask.any():
        unmatched_count = int(unmatched_feature_mask.sum())
        raise ValueError(
            f"{unmatched_count} feature rows could not be matched "
            "to the cell table."
        )

    merged = merged.drop(columns="_merge")

    group_columns = [sample_col, cell_type_col]

    summary = _aggregate_metrics(
        merged,
        group_columns=group_columns,
        metric_columns=selected_metrics,
    )

    cell_counts = (
        merged.groupby(
            group_columns,
            sort=False,
            dropna=False,
        )
        .size()
        .rename("cell_count")
        .reset_index()
    )

    sample_counts = (
        merged.groupby(
            sample_col,
            sort=False,
            dropna=False,
        )
        .size()
        .rename("sample_cell_count")
        .reset_index()
    )

    output = cell_counts.merge(
        sample_counts,
        on=sample_col,
        how="left",
        validate="many_to_one",
    )

    output["cell_fraction"] = (
        output["cell_count"] / output["sample_cell_count"]
    )

    output = output.merge(
        summary,
        on=group_columns,
        how="left",
        validate="one_to_one",
    )

    leading_columns = [
        sample_col,
        cell_type_col,
        "cell_count",
        "sample_cell_count",
        "cell_fraction",
    ]

    remaining_columns = [
        column
        for column in output.columns
        if column not in leading_columns
    ]

    return output[leading_columns + remaining_columns]


def aggregate_metrics_by_sample(
    features: pd.DataFrame,
    *,
    metric_columns: Sequence[str] | None = None,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
) -> pd.DataFrame:
    """
    Aggregate cell-level network metrics for each biological sample.

    Parameters
    ----------
    features:
        Cell-level output from ``compute_network_features``.
    metric_columns:
        Network metrics to aggregate. If omitted, all available default
        Layer 2 metrics are used.
    cell_id_col:
        Column containing cell identifiers.
    sample_col:
        Column containing sample identifiers.

    Returns
    -------
    pandas.DataFrame
        One row per sample containing cell count and descriptive statistics
        for every selected network metric.
    """

    working_features, selected_metrics = _prepare_features(
        features,
        metric_columns=metric_columns,
        cell_id_col=cell_id_col,
        sample_col=sample_col,
    )

    summary = _aggregate_metrics(
        working_features,
        group_columns=[sample_col],
        metric_columns=selected_metrics,
    )

    cell_counts = (
        working_features.groupby(
            sample_col,
            sort=False,
            dropna=False,
        )
        .size()
        .rename("cell_count")
        .reset_index()
    )

    output = cell_counts.merge(
        summary,
        on=sample_col,
        how="left",
        validate="one_to_one",
    )

    leading_columns = [
        sample_col,
        "cell_count",
    ]

    remaining_columns = [
        column
        for column in output.columns
        if column not in leading_columns
    ]

    return output[leading_columns + remaining_columns]
