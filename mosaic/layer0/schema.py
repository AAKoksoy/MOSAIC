"""
Core schema definitions for MOSAIC Layer 0.

This module defines the standard column names used by the MOSAIC
reference implementation. Users may still supply alternative column
names through function parameters.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CellTableSchema:
    """
    Standard column names for a MOSAIC cell table.

    Attributes
    ----------
    cell_id:
        Unique cell identifier within a sample.
    x:
        Cell x coordinate.
    y:
        Cell y coordinate.
    cell_type:
        Cell phenotype or cell-type label.
    sample_id:
        Sample or slide identifier.
    """

    cell_id: str = "cell_id"
    x: str = "x"
    y: str = "y"
    cell_type: str = "cell_type"
    sample_id: str = "sample_id"

    @property
    def required_columns(self) -> tuple[str, ...]:
        """
        Return the required Layer 0 column names.
        """

        return (
            self.cell_id,
            self.x,
            self.y,
            self.cell_type,
            self.sample_id,
        )


DEFAULT_CELL_SCHEMA = CellTableSchema()


OPTIONAL_CELL_COLUMNS = (
    "patient_id",
    "disease_type",
    "tissue_compartment",
    "imaging_platform",
    "clinical_outcome",
    "experimental_batch",
    "time_point",
)


def get_required_cell_columns(
    schema: CellTableSchema = DEFAULT_CELL_SCHEMA,
) -> tuple[str, ...]:
    """
    Return the required columns for a cell-table schema.

    Parameters
    ----------
    schema:
        Cell-table schema definition.

    Returns
    -------
    tuple[str, ...]
        Required column names.
    """

    return schema.required_columns
