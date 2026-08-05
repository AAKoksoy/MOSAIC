"""Regression tests for MOSAIC Layer 2 centrality features."""

import pandas as pd
import pytest

from mosaic.layer2.centrality import (
    compute_betweenness,
    compute_clustering_coefficient,
    compute_degree,
    compute_network_features,
    compute_pagerank,
    compute_weighted_degree,
)


def _path_graph_tables():
    """Return the path A-B-C plus isolated cell D."""

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
            "weight": [0.5, 0.25],
        }
    )

    return cells, edges


def test_compute_degree_for_path_and_isolated_cell():
    """A-B-C should have degrees 1, 2, 1; D should have degree 0."""

    cells, edges = _path_graph_tables()

    degree = compute_degree(cells, edges)
    values = degree.set_index("cell_id")["degree"].to_dict()

    assert values == {
        "A": 1,
        "B": 2,
        "C": 1,
        "D": 0,
    }


def test_compute_weighted_degree_sums_incident_weights():
    """Weighted degree should equal the sum of incident edge weights."""

    cells, edges = _path_graph_tables()

    weighted_degree = compute_weighted_degree(cells, edges)
    values = weighted_degree.set_index(
        "cell_id"
    )["weighted_degree"].to_dict()

    assert values["A"] == pytest.approx(0.5)
    assert values["B"] == pytest.approx(0.75)
    assert values["C"] == pytest.approx(0.25)
    assert values["D"] == pytest.approx(0.0)


def test_duplicate_edge_weights_are_combined():
    """Repeated undirected edges should contribute their combined weight."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "cell_id": ["A", "B"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "source": ["A", "B"],
            "target": ["B", "A"],
            "weight": [0.4, 0.6],
        }
    )

    weighted_degree = compute_weighted_degree(cells, edges)
    values = weighted_degree.set_index(
        "cell_id"
    )["weighted_degree"].to_dict()

    assert values["A"] == pytest.approx(1.0)
    assert values["B"] == pytest.approx(1.0)


def test_pagerank_is_calculated_independently_per_sample():
    """PageRank values should sum to one within each biological sample."""

    cells = pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_1",
                "sample_2",
                "sample_2",
            ],
            "cell_id": ["A", "B", "C", "A", "B"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": [
                "sample_1",
                "sample_1",
                "sample_2",
            ],
            "source": ["A", "B", "A"],
            "target": ["B", "C", "B"],
            "weight": [1.0, 1.0, 1.0],
        }
    )

    pagerank = compute_pagerank(cells, edges)

    sample_totals = pagerank.groupby(
        "sample_id"
    )["pagerank"].sum()

    assert sample_totals["sample_1"] == pytest.approx(1.0)
    assert sample_totals["sample_2"] == pytest.approx(1.0)

    sample_1 = pagerank.loc[
        pagerank["sample_id"] == "sample_1"
    ].set_index("cell_id")

    assert (
        sample_1.loc["B", "pagerank"]
        > sample_1.loc["A", "pagerank"]
    )
    assert sample_1.loc["A", "pagerank"] == pytest.approx(
        sample_1.loc["C", "pagerank"]
    )


def test_betweenness_identifies_middle_of_path():
    """B should lie on the only shortest path between A and C."""

    cells, edges = _path_graph_tables()

    betweenness = compute_betweenness(cells, edges)
    values = betweenness.set_index(
        "cell_id"
    )["betweenness"].to_dict()

    assert values["A"] == pytest.approx(0.0)
    assert values["B"] == pytest.approx(1 / 3)
    assert values["C"] == pytest.approx(0.0)
    assert values["D"] == pytest.approx(0.0)


def test_triangle_has_clustering_coefficient_one():
    """Every node in a complete three-node triangle should equal one."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 4,
            "cell_id": ["A", "B", "C", "D"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 3,
            "source": ["A", "B", "C"],
            "target": ["B", "C", "A"],
        }
    )

    clustering = compute_clustering_coefficient(cells, edges)
    values = clustering.set_index(
        "cell_id"
    )["clustering_coefficient"].to_dict()

    assert values["A"] == pytest.approx(1.0)
    assert values["B"] == pytest.approx(1.0)
    assert values["C"] == pytest.approx(1.0)
    assert values["D"] == pytest.approx(0.0)


def test_path_graph_has_zero_clustering():
    """A path contains no closed triangles."""

    cells, edges = _path_graph_tables()

    clustering = compute_clustering_coefficient(cells, edges)

    assert (
        clustering["clustering_coefficient"] == pytest.approx(0.0)
    ).all()


def test_compute_network_features_returns_complete_table():
    """The combined function should return every Layer 2 metric."""

    cells, edges = _path_graph_tables()

    features = compute_network_features(cells, edges)

    assert list(features.columns) == [
        "sample_id",
        "cell_id",
        "degree",
        "weighted_degree",
        "pagerank",
        "betweenness",
        "clustering_coefficient",
    ]

    assert len(features) == len(cells)
    assert not features.isna().any().any()

    isolated = features.loc[
        features["cell_id"] == "D"
    ].iloc[0]

    assert isolated["degree"] == 0
    assert isolated["weighted_degree"] == pytest.approx(0.0)
    assert isolated["betweenness"] == pytest.approx(0.0)
    assert isolated["clustering_coefficient"] == pytest.approx(0.0)


def test_repeated_cell_ids_remain_separate_across_samples():
    """Identical IDs in different samples must receive separate metrics."""

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
        }
    )

    degree = compute_degree(cells, edges)

    sample_1 = degree.loc[
        degree["sample_id"] == "sample_1"
    ].set_index("cell_id")["degree"]

    sample_2 = degree.loc[
        degree["sample_id"] == "sample_2"
    ].set_index("cell_id")["degree"]

    assert sample_1.to_dict() == {"A": 1, "B": 1}
    assert sample_2.to_dict() == {"A": 0, "B": 0}
