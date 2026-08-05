"""Regression tests for MOSAIC Layer 2 connected components."""

import pandas as pd
import pytest

from mosaic.layer2.components import (
    compute_connected_components,
    connected_component_summary,
)


def test_compute_connected_components_preserves_isolated_cells():
    """A-B-C should form one component, while D remains isolated."""

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
        }
    )

    components = compute_connected_components(cells, edges)

    connected = components.loc[
        components["cell_id"].isin(["A", "B", "C"])
    ]
    isolated = components.loc[
        components["cell_id"] == "D"
    ].iloc[0]

    assert len(components) == 4
    assert connected["component_id"].nunique() == 1
    assert set(connected["component_size"]) == {3}
    assert not connected["is_isolated_component"].any()

    assert isolated["component_size"] == 1
    assert bool(isolated["is_isolated_component"])


def test_component_ids_are_deterministic():
    """Components should receive stable IDs based on size and node name."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 5,
            "cell_id": ["E", "D", "C", "B", "A"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "source": ["B", "D"],
            "target": ["A", "E"],
        }
    )

    first_result = compute_connected_components(cells, edges)
    second_result = compute_connected_components(
        cells.sample(frac=1.0, random_state=42),
        edges.sample(frac=1.0, random_state=42),
    )

    first_mapping = (
        first_result.set_index("cell_id")["component_id"].to_dict()
    )
    second_mapping = (
        second_result.set_index("cell_id")["component_id"].to_dict()
    )

    assert first_mapping == second_mapping
    assert first_mapping["A"] == "sample_1_component_1"
    assert first_mapping["B"] == "sample_1_component_1"
    assert first_mapping["D"] == "sample_1_component_2"
    assert first_mapping["E"] == "sample_1_component_2"
    assert first_mapping["C"] == "sample_1_component_3"


def test_unknown_edge_nodes_are_excluded():
    """Edges to cells absent from the cell table must not create nodes."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "cell_id": ["A", "B"],
        }
    )

    edges = pd.DataFrame(
        {
            "sample_id": ["sample_1"],
            "source": ["A"],
            "target": ["UNKNOWN"],
        }
    )

    components = compute_connected_components(cells, edges)

    assert set(components["cell_id"]) == {"A", "B"}
    assert len(components) == 2
    assert set(components["component_size"]) == {1}
    assert components["is_isolated_component"].all()


def test_connected_component_summary_returns_expected_statistics():
    """Component summaries should reflect sizes three, one, and one."""

    components = pd.DataFrame(
        {
            "sample_id": ["sample_1"] * 5,
            "cell_id": ["A", "B", "C", "D", "E"],
            "component_id": [
                "sample_1_component_1",
                "sample_1_component_1",
                "sample_1_component_1",
                "sample_1_component_2",
                "sample_1_component_3",
            ],
            "component_size": [3, 3, 3, 1, 1],
            "is_isolated_component": [
                False,
                False,
                False,
                True,
                True,
            ],
        }
    )

    summary = connected_component_summary(components)

    assert len(summary) == 1

    result = summary.iloc[0]

    assert result["sample_id"] == "sample_1"
    assert result["component_count"] == 3
    assert result["largest_component_size"] == 3
    assert result["largest_component_fraction"] == pytest.approx(3 / 5)
    assert result["isolated_component_count"] == 2
    assert result["mean_component_size"] == pytest.approx(5 / 3)
    assert result["median_component_size"] == pytest.approx(1.0)


def test_repeated_cell_ids_remain_separate_across_samples():
    """Identical cell IDs in different samples must not be connected."""

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

    components = compute_connected_components(cells, edges)

    sample_1 = components.loc[
        components["sample_id"] == "sample_1"
    ]
    sample_2 = components.loc[
        components["sample_id"] == "sample_2"
    ]

    assert sample_1["component_id"].nunique() == 1
    assert set(sample_1["component_size"]) == {2}
    assert not sample_1["is_isolated_component"].any()

    assert sample_2["component_id"].nunique() == 2
    assert set(sample_2["component_size"]) == {1}
    assert sample_2["is_isolated_component"].all()


def test_empty_component_table_returns_expected_columns():
    """An empty component table should produce an empty valid summary."""

    components = pd.DataFrame(
        columns=[
            "sample_id",
            "component_id",
            "component_size",
        ]
    )

    summary = connected_component_summary(components)

    assert summary.empty
    assert list(summary.columns) == [
        "sample_id",
        "component_count",
        "largest_component_size",
        "largest_component_fraction",
        "isolated_component_count",
        "mean_component_size",
        "median_component_size",
    ]
