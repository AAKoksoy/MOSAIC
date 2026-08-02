"""
Connected-component analysis for MOSAIC Layer 2.

Components are calculated independently for each biological sample to prevent
cells with identical identifiers in different samples from being mixed.
"""

from __future__ import annotations

import networkx as nx
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


def compute_connected_components(
    cells: pd.DataFrame,
    edges: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    sample_col: str = "sample_id",
    source_col: str = "source",
    target_col: str = "target",
) -> pd.DataFrame:
    """
    Assign every cell to a connected component within its sample.

    All cells are added to the graph before edges are added. Therefore,
    isolated cells are retained as one-node components.

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
        One row per cell with:

        - sample identifier,
        - cell identifier,
        - component identifier,
        - component size,
        - isolated-component indicator.
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

    if cells.empty:
        raise ValueError("cells cannot be empty.")

    if cells[[cell_id_col, sample_col]].isna().any().any():
        raise ValueError(
            "Cell identifiers and sample identifiers cannot be missing."
        )

    working_cells = cells[[sample_col, cell_id_col]].copy()
    working_cells[sample_col] = working_cells[sample_col].astype(str)
    working_cells[cell_id_col] = working_cells[cell_id_col].astype(str)

    if working_cells.duplicated(
        subset=[sample_col, cell_id_col]
    ).any():
        raise ValueError(
            "Cell identifiers must be unique within each sample."
        )

    working_edges = edges[
        [sample_col, source_col, target_col]
    ].copy()

    working_edges[sample_col] = (
        working_edges[sample_col].astype(str)
    )
    working_edges[source_col] = (
        working_edges[source_col].astype(str)
    )
    working_edges[target_col] = (
        working_edges[target_col].astype(str)
    )

    component_records: list[dict[str, object]] = []

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

        valid_edges: list[tuple[str, str]] = []

        for source_value, target_value in sample_edges[
            [source_col, target_col]
        ].itertuples(index=False, name=None):
            if (
                source_value in known_nodes
                and target_value in known_nodes
            ):
                valid_edges.append(
                    (source_value, target_value)
                )

        graph.add_edges_from(valid_edges)

        components = list(nx.connected_components(graph))

        components.sort(
            key=lambda component: (
                -len(component),
                min(str(node) for node in component),
            )
        )

        for component_number, component_nodes in enumerate(
            components,
            start=1,
        ):
            component_id = (
                f"{sample_value}_component_{component_number}"
            )
            component_size = len(component_nodes)

            for node in sorted(
                component_nodes,
                key=str,
            ):
                component_records.append(
                    {
                        sample_col: sample_value,
                        cell_id_col: node,
                        "component_id": component_id,
                        "component_size": component_size,
                        "is_isolated_component": (
                            component_size == 1
                        ),
                    }
                )

    output = pd.DataFrame(component_records)

    output_columns = [
        sample_col,
        cell_id_col,
        "component_id",
        "component_size",
        "is_isolated_component",
    ]

    return output[output_columns]


def connected_component_summary(
    components: pd.DataFrame,
    *,
    sample_col: str = "sample_id",
    component_col: str = "component_id",
    component_size_col: str = "component_size",
) -> pd.DataFrame:
    """
    Summarize connected-component organization for each sample.

    Parameters
    ----------
    components:
        Output from ``compute_connected_components``.
    sample_col:
        Column containing sample identifiers.
    component_col:
        Column containing component identifiers.
    component_size_col:
        Column containing component sizes.

    Returns
    -------
    pandas.DataFrame
        One row per sample containing component counts, largest-component
        statistics, isolated-component counts, and size summaries.
    """

    _require_dataframe(components, "components")

    _require_columns(
        components,
        [
            sample_col,
            component_col,
            component_size_col,
        ],
        "components",
    )

    output_columns = [
        sample_col,
        "component_count",
        "largest_component_size",
        "largest_component_fraction",
        "isolated_component_count",
        "mean_component_size",
        "median_component_size",
    ]

    if components.empty:
        return pd.DataFrame(columns=output_columns)

    working_components = components.copy()
    working_components[sample_col] = (
        working_components[sample_col].astype(str)
    )

    component_table = (
        working_components[
            [
                sample_col,
                component_col,
                component_size_col,
            ]
        ]
        .drop_duplicates(
            subset=[sample_col, component_col]
        )
        .copy()
    )

    summary_records: list[dict[str, object]] = []

    for sample_value, sample_components in component_table.groupby(
        sample_col,
        sort=False,
        dropna=False,
    ):
        component_sizes = pd.to_numeric(
            sample_components[component_size_col],
            errors="raise",
        )

        if (component_sizes < 1).any():
            raise ValueError(
                "Component sizes must be greater than or equal to one."
            )

        node_count = int(component_sizes.sum())
        largest_component_size = int(component_sizes.max())

        summary_records.append(
            {
                sample_col: sample_value,
                "component_count": int(len(sample_components)),
                "largest_component_size": largest_component_size,
                "largest_component_fraction": float(
                    largest_component_size / node_count
                ),
                "isolated_component_count": int(
                    (component_sizes == 1).sum()
                ),
                "mean_component_size": float(
                    component_sizes.mean()
                ),
                "median_component_size": float(
                    component_sizes.median()
                ),
            }
        )

    return pd.DataFrame(
        summary_records,
        columns=output_columns,
    )
