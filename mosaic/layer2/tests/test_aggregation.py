"""Regression tests for MOSAIC Layer 2 feature aggregation."""

import pandas as pd
import pytest

from mosaic.layer2.aggregation import (
    aggregate_metrics_by_cell_type,
    aggregate_metrics_by_sample,
)


def _feature_table():
    """Return known cell-level metrics for two samples."""

    return pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_2",
                "sample_2",
            ],
            "cell_id": ["A", "B", "C", "D", "A", "B"],
            "degree": [1, 3, 2, 0, 1, 1],
            "weighted_degree": [0.5, 1.5, 1.0, 0.0, 0.8, 0.8],
            "pagerank": [0.20, 0.40, 0.30, 0.10, 0.50, 0.50],
            "betweenness": [0.0, 0.5, 0.2, 0.0, 0.0, 0.0],
            "clustering_coefficient": [
                0.0,
                0.5,
                1.0,
                0.0,
                0.0,
                0.0,
            ],
        }
    )


def _cell_table():
    """Return cell-type annotations corresponding to the feature table."""

    return pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_2",
                "sample_2",
            ],
            "cell_id": ["A", "B", "C", "D", "A", "B"],
            "cell_type": [
                "T_cell",
                "T_cell",
                "Tumor",
                "Tumor",
                "T_cell",
                "Tumor",
            ],
        }
    )


def test_aggregate_metrics_by_sample_returns_expected_counts():
    """Each sample should have the correct number of cells."""

    summary = aggregate_metrics_by_sample(_feature_table())

    counts = summary.set_index("sample_id")["cell_count"].to_dict()

    assert counts == {
        "sample_1": 4,
        "sample_2": 2,
    }


def test_sample_aggregation_returns_expected_statistics():
    """Known degree values should produce known summary statistics."""

    summary = aggregate_metrics_by_sample(_feature_table())

    sample_1 = summary.loc[
        summary["sample_id"] == "sample_1"
    ].iloc[0]

    assert sample_1["degree_mean"] == pytest.approx(1.5)
    assert sample_1["degree_median"] == pytest.approx(1.5)
    assert sample_1["degree_std"] == pytest.approx(
        1.2909944487358056
    )
    assert sample_1["degree_min"] == 0
    assert sample_1["degree_max"] == 3

    assert sample_1["weighted_degree_mean"] == pytest.approx(0.75)
    assert sample_1["pagerank_mean"] == pytest.approx(0.25)
    assert sample_1["betweenness_max"] == pytest.approx(0.5)
    assert (
        sample_1["clustering_coefficient_max"]
        == pytest.approx(1.0)
    )


def test_cell_type_aggregation_returns_counts_and_fractions():
    """Cell-type counts and within-sample fractions should be correct."""

    summary = aggregate_metrics_by_cell_type(
        _feature_table(),
        _cell_table(),
    )

    sample_1_t_cells = summary.loc[
        (summary["sample_id"] == "sample_1")
        & (summary["cell_type"] == "T_cell")
    ].iloc[0]

    assert sample_1_t_cells["cell_count"] == 2
    assert sample_1_t_cells["sample_cell_count"] == 4
    assert sample_1_t_cells["cell_fraction"] == pytest.approx(0.5)

    sample_2_tumor = summary.loc[
        (summary["sample_id"] == "sample_2")
        & (summary["cell_type"] == "Tumor")
    ].iloc[0]

    assert sample_2_tumor["cell_count"] == 1
    assert sample_2_tumor["sample_cell_count"] == 2
    assert sample_2_tumor["cell_fraction"] == pytest.approx(0.5)


def test_cell_type_aggregation_returns_expected_metrics():
    """Metrics should be summarized within sample and cell type."""

    summary = aggregate_metrics_by_cell_type(
        _feature_table(),
        _cell_table(),
    )

    sample_1_t_cells = summary.loc[
        (summary["sample_id"] == "sample_1")
        & (summary["cell_type"] == "T_cell")
    ].iloc[0]

    assert sample_1_t_cells["degree_mean"] == pytest.approx(2.0)
    assert sample_1_t_cells["degree_median"] == pytest.approx(2.0)
    assert sample_1_t_cells["degree_std"] == pytest.approx(
        1.4142135623730951
    )
    assert sample_1_t_cells["degree_min"] == 1
    assert sample_1_t_cells["degree_max"] == 3

    assert (
        sample_1_t_cells["weighted_degree_mean"]
        == pytest.approx(1.0)
    )
    assert sample_1_t_cells["pagerank_mean"] == pytest.approx(0.30)
    assert sample_1_t_cells["betweenness_mean"] == pytest.approx(
        0.25
    )
    assert (
        sample_1_t_cells["clustering_coefficient_mean"]
        == pytest.approx(0.25)
    )


def test_single_cell_group_has_zero_standard_deviation():
    """A one-cell group should report zero rather than missing std."""

    summary = aggregate_metrics_by_cell_type(
        _feature_table(),
        _cell_table(),
    )

    sample_2_t_cells = summary.loc[
        (summary["sample_id"] == "sample_2")
        & (summary["cell_type"] == "T_cell")
    ].iloc[0]

    standard_deviation_columns = [
        column
        for column in summary.columns
        if column.endswith("_std")
    ]

    assert (
        sample_2_t_cells[standard_deviation_columns] == 0.0
    ).all()


def test_repeated_cell_ids_are_joined_using_sample_and_cell_id():
    """Repeated cell IDs must receive the annotation from their sample."""

    features = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "cell_id": ["A", "A"],
            "degree": [1.0, 5.0],
        }
    )

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "cell_id": ["A", "A"],
            "cell_type": ["T_cell", "Tumor"],
        }
    )

    summary = aggregate_metrics_by_cell_type(
        features,
        cells,
        metric_columns=["degree"],
    )

    sample_1 = summary.loc[
        summary["sample_id"] == "sample_1"
    ].iloc[0]

    sample_2 = summary.loc[
        summary["sample_id"] == "sample_2"
    ].iloc[0]

    assert sample_1["cell_type"] == "T_cell"
    assert sample_1["degree_mean"] == pytest.approx(1.0)

    assert sample_2["cell_type"] == "Tumor"
    assert sample_2["degree_mean"] == pytest.approx(5.0)


def test_selected_metric_columns_limit_output():
    """Only explicitly selected metrics should be aggregated."""

    summary = aggregate_metrics_by_sample(
        _feature_table(),
        metric_columns=["degree", "pagerank"],
    )

    assert "degree_mean" in summary.columns
    assert "pagerank_mean" in summary.columns
    assert "weighted_degree_mean" not in summary.columns
    assert "betweenness_mean" not in summary.columns
    assert "clustering_coefficient_mean" not in summary.columns


def test_cell_type_fractions_sum_to_one_within_each_sample():
    """All represented cell-type fractions should total one per sample."""

    summary = aggregate_metrics_by_cell_type(
        _feature_table(),
        _cell_table(),
    )

    fraction_totals = summary.groupby(
        "sample_id"
    )["cell_fraction"].sum()

    assert fraction_totals["sample_1"] == pytest.approx(1.0)
    assert fraction_totals["sample_2"] == pytest.approx(1.0)
