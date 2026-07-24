"""
Layer 1 — Spatial Graph Construction

This module converts validated cell-level observations into spatial graphs.
Cells are represented as nodes, and spatial relationships are represented
as edges according to explicit neighborhood rules.
"""

from .radius_graph import (
    build_radius_graph,
    summarize_spatial_graph,
    validate_edge_table,
)

__all__ = [
    "build_radius_graph",
    "validate_edge_table",
    "summarize_spatial_graph",
]
