"""Regression tests for MOSAIC Layer 2 input validation."""

import numpy as np
import pandas as pd
import pytest

from mosaic.layer2.aggregation import (
    aggregate_metrics_by_cell_type,
    aggregate_metrics_by_sample,
)
from mosaic.layer2.centrality import (
    compute_betweenness,
    compute_degree,
    compute_pagerank,
    compute_weighted_degree,
)
from mosaic.layer2.components import compute_connected_components
from mosaic.layer2.qc import calculate_graph_qc


def _valid_cells():
    """Return a minimal valid cell table."""

    return pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "cell_id": ["A", "B"],
        }
    )


def _valid_edges():
    """Return a minimal valid edge table for all Layer 2 functions."""

    return pd.DataFrame(
        {
            "sample_id": ["sample_1"],
            "source": ["A"],
            "target": ["B"],
            "distance_um": [10.0],
            "weight": [0.5],
        }
    )


def test_missing_required_cell_column_raises_clear_error():
    """Centrality functions should reject incomplete cell tables."""

    cells = pd.DataFrame({"cell_id": ["A", "B"]})

    with pytest.raises(
        ValueError,
        match="Missing required columns in cells: sample_id",
    ):
        compute_degree(cells, _valid_edges())


def test_qc_requires_distance_column():
    """Graph QC should require spatial edge distances."""

    edges = _valid_edges().drop(columns="distance_um")

    with pytest.raises(
        ValueError,
        match="Missing required columns in edges: distance_um",
    ):
        calculate_graph_qc(_valid_cells(), edges)


def test_empty_cell_table_is_rejected():
    """A graph cannot be calculated without biological entities."""

    empty_cells = pd.DataFrame(
        columns=["sample_id", "cell_id"]
    )

    with pytest.raises(
        ValueError,
        match="cells cannot be empty",
    ):
        compute_connected_components(
            empty_cells,
            _valid_edges(),
        )


def test_duplicate_cell_ids_within_sample_are_rejected():
    """Cell IDs must be unique inside each biological sample."""

    cells = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "cell_id": ["A", "A"],
        }
    )

    with pytest.raises(
        ValueError,
        match=(
            "Cell identifiers must be unique "
            "within each sample"
        ),
    ):
        compute_degree(cells, _valid_edges())


@pytest.mark.parametrize(
    "invalid_weight",
    [-1.0, np.nan, np.inf, "bad"],
)
def test_invalid_edge_weights_are_rejected(invalid_weight):
    """Weights must be finite, non-negative numeric values."""

    edges = _valid_edges()
    edges["weight"] = [invalid_weight]

    with pytest.raises(
        ValueError,
        match=(
            "weight must contain finite, "
            "non-negative numeric values"
        ),
    ):
        compute_weighted_degree(
            _valid_cells(),
            edges,
        )


@pytest.mark.parametrize(
    "invalid_alpha",
    [0.0, 1.0, -0.1, np.inf],
)
def test_invalid_pagerank_alpha_values_are_rejected(
    invalid_alpha,
):
    """PageRank alpha must be finite and strictly between 0 and 1."""

    with pytest.raises(
        ValueError,
        match=(
            "alpha must be greater than zero "
            "and less than one"
        ),
    ):
        compute_pagerank(
            _valid_cells(),
            _valid_edges(),
            alpha=invalid_alpha,
        )


def test_non_numeric_pagerank_alpha_is_rejected():
    """PageRank should distinguish bad types from bad ranges."""

    with pytest.raises(
        TypeError,
        match="alpha must be numeric",
    ):
        compute_pagerank(
            _valid_cells(),
            _valid_edges(),
            alpha="0.85",
        )


@pytest.mark.parametrize(
    "invalid_k",
    [0, -1],
)
def test_non_positive_approximate_betweenness_k_is_rejected(
    invalid_k,
):
    """Approximate betweenness must sample at least one node."""

    with pytest.raises(
        ValueError,
        match=(
            "approximate_k must be greater than "
            "or equal to one"
        ),
    ):
        compute_betweenness(
            _valid_cells(),
            _valid_edges(),
            approximate_k=invalid_k,
        )


@pytest.mark.parametrize(
    "invalid_k",
    [1.5, "2", True],
)
def test_invalid_approximate_betweenness_k_types_are_rejected(
    invalid_k,
):
    """Approximate sample size must be an integer or None."""

    with pytest.raises(
        TypeError,
        match=(
            "approximate_k must be an integer or None"
        ),
    ):
        compute_betweenness(
            _valid_cells(),
            _valid_edges(),
            approximate_k=invalid_k,
        )


def test_unmatched_feature_rows_are_rejected():
    """Every feature must match a sample-aware cell annotation."""

    features = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_1"],
            "cell_id": ["A", "UNKNOWN"],
            "degree": [1.0, 2.0],
        }
    )

    cells = _valid_cells().assign(
        cell_type=["T_cell", "Tumor"]
    )

    with pytest.raises(
        ValueError,
        match=(
            "1 feature rows could not be matched "
            "to the cell table"
        ),
    ):
        aggregate_metrics_by_cell_type(
            features,
            cells,
            metric_columns=["degree"],
        )


def test_metric_columns_cannot_be_a_single_string():
    """Metric selection must be a sequence, not a bare string."""

    features = pd.DataFrame(
        {
            "sample_id": ["sample_1"],
            "cell_id": ["A"],
            "degree": [1.0],
        }
    )

    with pytest.raises(
        TypeError,
        match=(
            "metric_columns must be a sequence "
            "of column names"
        ),
    ):
        aggregate_metrics_by_sample(
            features,
            metric_columns="degree",
        )


@pytest.mark.parametrize(
    "invalid_value",
    [np.nan, np.inf, "bad"],
)
def test_invalid_network_metric_values_are_rejected(
    invalid_value,
):
    """Aggregation should reject invalid feature values."""

    features = pd.DataFrame(
        {
            "sample_id": ["sample_1"],
            "cell_id": ["A"],
            "degree": [invalid_value],
        }
    )

    with pytest.raises(
        ValueError,
        match="degree must contain finite numeric values",
    ):
        aggregate_metrics_by_sample(
            features,
            metric_columns=["degree"],
        )
