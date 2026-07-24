"""
Validation utilities for MOSAIC Layer 0.

Layer 0 standardizes and validates cell-level observations before
spatial graph construction.
"""

from __future__ import annotations

from typing import Any

import pandas as pd


DEFAULT_REQUIRED_COLUMNS = (
    "cell_id",
    "x",
    "y",
    "cell_type",
    "sample_id",
)


def validate_cell_table(
    cells: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    x_col: str = "x",
    y_col: str = "y",
    phenotype_col: str = "cell_type",
    sample_col: str = "sample_id",
    allow_missing_phenotypes: bool = False,
) -> pd.DataFrame:
    """
    Validate a cell-level spatial biology table.

    Parameters
    ----------
    cells:
        Input cell table containing one row per cell.
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
    allow_missing_phenotypes:
        Whether missing phenotype values are permitted.

    Returns
    -------
    pandas.DataFrame
        A validated copy of the input table.

    Raises
    ------
    TypeError
        If the input is not a pandas DataFrame.
    ValueError
        If required columns are missing or invalid values are detected.

    Notes
    -----
    Cell identifiers are required to be unique within each sample.
    The same cell identifier may appear in different samples.
    """

    if not isinstance(cells, pd.DataFrame):
        raise TypeError("cells must be a pandas DataFrame.")

    required_columns = [
        cell_id_col,
        x_col,
        y_col,
        phenotype_col,
        sample_col,
    ]

    missing_columns = [
        column for column in required_columns if column not in cells.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(str(column) for column in missing_columns)
        )

    if cells.empty:
        raise ValueError("The cell table is empty.")

    validated = cells.copy()

    identifier_columns = [cell_id_col, sample_col]

    for column in identifier_columns:
        missing_count = int(validated[column].isna().sum())

        if missing_count > 0:
            raise ValueError(
                f"Column '{column}' contains {missing_count} missing value(s)."
            )

    coordinate_columns = [x_col, y_col]

    for column in coordinate_columns:
        converted = pd.to_numeric(validated[column], errors="coerce")

        invalid_mask = converted.isna()

        if invalid_mask.any():
            invalid_count = int(invalid_mask.sum())

            raise ValueError(
                f"Column '{column}' contains {invalid_count} missing or "
                "nonnumeric value(s)."
            )

        validated[column] = converted.astype(float)

    if not allow_missing_phenotypes:
        missing_phenotypes = int(validated[phenotype_col].isna().sum())

        if missing_phenotypes > 0:
            raise ValueError(
                f"Column '{phenotype_col}' contains "
                f"{missing_phenotypes} missing value(s)."
            )

    duplicate_mask = validated.duplicated(
        subset=[sample_col, cell_id_col],
        keep=False,
    )

    if duplicate_mask.any():
        duplicate_count = int(duplicate_mask.sum())

        raise ValueError(
            f"Found {duplicate_count} row(s) with duplicate cell identifiers "
            "within the same sample."
        )

    return validated


def standardize_cell_table(
    cells: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    x_col: str = "x",
    y_col: str = "y",
    phenotype_col: str = "cell_type",
    sample_col: str = "sample_id",
    strip_text: bool = True,
) -> pd.DataFrame:
    """
    Standardize core columns in a cell-level table.

    The function preserves all original columns while standardizing
    identifiers, coordinates, phenotype labels, and sample labels.

    Parameters
    ----------
    cells:
        Input cell table.
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
    strip_text:
        Whether surrounding whitespace should be removed from text fields.

    Returns
    -------
    pandas.DataFrame
        A standardized copy of the cell table.
    """

    standardized = validate_cell_table(
        cells,
        cell_id_col=cell_id_col,
        x_col=x_col,
        y_col=y_col,
        phenotype_col=phenotype_col,
        sample_col=sample_col,
    )

    standardized[cell_id_col] = standardized[cell_id_col].astype(str)
    standardized[sample_col] = standardized[sample_col].astype(str)
    standardized[phenotype_col] = standardized[phenotype_col].astype(str)

    if strip_text:
        standardized[cell_id_col] = standardized[cell_id_col].str.strip()
        standardized[sample_col] = standardized[sample_col].str.strip()
        standardized[phenotype_col] = standardized[phenotype_col].str.strip()

    empty_cell_ids = standardized[cell_id_col].eq("")
    empty_sample_ids = standardized[sample_col].eq("")
    empty_phenotypes = standardized[phenotype_col].eq("")

    if empty_cell_ids.any():
        raise ValueError("Cell identifiers cannot be empty strings.")

    if empty_sample_ids.any():
        raise ValueError("Sample identifiers cannot be empty strings.")

    if empty_phenotypes.any():
        raise ValueError("Phenotype labels cannot be empty strings.")

    standardized[x_col] = standardized[x_col].astype(float)
    standardized[y_col] = standardized[y_col].astype(float)

    return standardized


def report_cell_qc(
    cells: pd.DataFrame,
    *,
    cell_id_col: str = "cell_id",
    x_col: str = "x",
    y_col: str = "y",
    phenotype_col: str = "cell_type",
    sample_col: str = "sample_id",
) -> dict[str, Any]:
    """
    Generate a quality-control summary for a cell table.

    Parameters
    ----------
    cells:
        Input cell table.
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

    Returns
    -------
    dict
        Summary statistics describing the cell table.
    """

    if not isinstance(cells, pd.DataFrame):
        raise TypeError("cells must be a pandas DataFrame.")

    required_columns = [
        cell_id_col,
        x_col,
        y_col,
        phenotype_col,
        sample_col,
    ]

    missing_columns = [
        column for column in required_columns if column not in cells.columns
    ]

    if missing_columns:
        raise ValueError(
            "Missing required columns: "
            + ", ".join(str(column) for column in missing_columns)
        )

    duplicate_count = int(
        cells.duplicated(
            subset=[sample_col, cell_id_col],
            keep=False,
        ).sum()
    )

    numeric_x = pd.to_numeric(cells[x_col], errors="coerce")
    numeric_y = pd.to_numeric(cells[y_col], errors="coerce")

    report = {
        "number_of_rows": int(len(cells)),
        "number_of_samples": int(cells[sample_col].nunique(dropna=True)),
        "number_of_cell_types": int(
            cells[phenotype_col].nunique(dropna=True)
        ),
        "missing_cell_ids": int(cells[cell_id_col].isna().sum()),
        "missing_sample_ids": int(cells[sample_col].isna().sum()),
        "missing_phenotypes": int(cells[phenotype_col].isna().sum()),
        "missing_or_invalid_x_coordinates": int(numeric_x.isna().sum()),
        "missing_or_invalid_y_coordinates": int(numeric_y.isna().sum()),
        "duplicate_cell_rows_within_sample": duplicate_count,
        "minimum_x": _safe_min(numeric_x),
        "maximum_x": _safe_max(numeric_x),
        "minimum_y": _safe_min(numeric_y),
        "maximum_y": _safe_max(numeric_y),
    }

    return report


def _safe_min(values: pd.Series) -> float | None:
    """Return the minimum numeric value or None when unavailable."""

    valid_values = values.dropna()

    if valid_values.empty:
        return None

    return float(valid_values.min())


def _safe_max(values: pd.Series) -> float | None:
    """Return the maximum numeric value or None when unavailable."""

    valid_values = values.dropna()

    if valid_values.empty:
        return None

    return float(valid_values.max())
