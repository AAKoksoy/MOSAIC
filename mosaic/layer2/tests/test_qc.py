"""Regression tests for MOSAIC Layer 2 graph quality control."""

import numpy as np
import pandas as pd
import pytest

from mosaic.layer2.qc import (
    calculate_graph_qc,
    find_high_degree_outliers,
    find_isolated_nodes,
)


def test_find_isolated_nodes_retains_sample_identity():
    """A cell with no incident edges should be reported as isolated."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 4,
            "cell_id": ["A", "B", "C", "D"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "source": ["A", "B"],
            "target": ["B", "C"],
            "distance_um": [10.0, 20.0],
        }
    )

    isolated = find_isolated_nodes(cells, edges)

    assert len(isolated) == 1
    assert isolated.iloc[0]["sample_id"] == "sample_1"
    assert isolated.iloc[0]["cell_id"] == "D"
    assert isolated.iloc[0]["degree"] == 0


def test_find_high_degree_outliers_detects_star_center():
    """The center of a sufficiently large star should be a degree outlier."""

    peripheral_nodes = [f"P{number}" for number in range(1, 9)]

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 9,
            "cell_id": ["CENTER"] + peripheral_nodes,
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 8,
            "source": ["CENTER"] * 8,
            "target": peripheral_nodes,
            "distance_um": [10.0] * 8,
        }
    )

    outliers = find_high_degree_outliers(cells, edges)

    assert len(outliers) == 1
    assert outliers.iloc[0]["sample_id"] == "sample_1"
    assert outliers.iloc[0]["cell_id"] == "CENTER"
    assert outliers.iloc[0]["degree"] == 8


def test_calculate_graph_qc_returns_expected_summary():
    """QC statistics should match a graph with one isolated cell."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 4,
            "cell_id": ["A", "B", "C", "D"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "source": ["A", "B"],
            "target": ["B", "C"],
            "distance_um": [10.0, 20.0],
        }
    )

    qc = calculate_graph_qc(cells, edges)

    assert len(qc) == 1

    result = qc.iloc[0]

    assert result["sample_id"] == "sample_1"
    assert result["node_count"] == 4
    assert result["edge_count"] == 2
    assert result["graph_density"] == pytest.approx(1 / 3)

    assert result["mean_degree"] == pytest.approx(1.0)
    assert result["median_degree"] == pytest.approx(1.0)
    assert result["maximum_degree"] == 2

    assert result["isolated_node_count"] == 1
    assert result["isolated_node_fraction"] == pytest.approx(0.25)
    assert result["high_degree_outlier_count"] == 0

    assert result["self_loop_count"] == 0
    assert result["duplicate_edge_count"] == 0
    assert result["missing_node_edge_count"] == 0
    assert result["cross_sample_violation_count"] == 0
    assert result["invalid_distance_count"] == 0

    assert result["minimum_distance_um"] == pytest.approx(10.0)
    assert result["mean_distance_um"] == pytest.approx(15.0)
    assert result["median_distance_um"] == pytest.approx(15.0)
    assert result["maximum_distance_um"] == pytest.approx(20.0)


def test_calculate_graph_qc_detects_edge_problems():
    """QC should identify malformed or biologically invalid edges."""

    cells = pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_2",
            ],
            "cell_id": ["A", "B", "A"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_1",
            ],
            "source": ["A", "A", "A", "A"],
            "target": ["A", "B", "B", "UNKNOWN"],
            "distance_um": [0.0, 10.0, np.nan, 5.0],
        }
    )

    qc = calculate_graph_qc(cells, edges)

    sample_1 = qc.loc[
        qc["sample_id"] == "sample_1"
    ].iloc[0]

    assert sample_1["self_loop_count"] == 1
    assert sample_1["duplicate_edge_count"] == 1
    assert sample_1["missing_node_edge_count"] == 1
    assert sample_1["invalid_distance_count"] == 1


def test_repeated_cell_ids_in_different_samples_remain_separate():
    """The same cell ID may occur safely in different samples."""

    cells = pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_2",
                "sample_2",
            ],
            "cell_id": ["A", "B", "A", "B"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1"],
            "source": ["A"],
            "target": ["B"],
            "distance_um": [10.0],
        }
    )

    isolated = find_isolated_nodes(cells, edges)

    sample_1_isolated = isolated.loc[
        isolated["sample_id"] == "sample_1"
    ]

    sample_2_isolated = isolated.loc[
        isolated["sample_id"] == "sample_2"
    ]

    assert sample_1_isolated.empty
    assert set(sample_2_isolated["cell_id"]) == {"A", "B"}
