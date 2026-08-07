# Changelog

All notable changes to MOSAIC Spatial will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/).


## [0.1.0] - 2026-08-07

### Added

- Layer 0: biological entity preparation and validation.
- Layer 1: sample-aware spatial graph construction.
- Radius-based neighborhood graph generation using spatial coordinates.
- Layer 2: cell-level network feature calculation.
- Degree and weighted-degree centrality.
- PageRank centrality.
- Betweenness centrality.
- Local clustering coefficient.
- Connected-component detection and summaries.
- Graph quality-control metrics.
- Isolated-cell and degree-outlier detection.
- Sample-level network feature aggregation.
- Cell-type-level network feature aggregation.
- Protection against incorrectly merging repeated cell identifiers across samples.
- Input validation with clear errors for malformed tables and invalid parameters.
- Regression tests covering MOSAIC Spatial Layers 0–2.
- Automated testing on Python 3.10, 3.11, and 3.12 through GitHub Actions.
- Python packaging and dependency configuration through `pyproject.toml`.

### Validation

- 52 regression tests passing.
- Complete test suite verified on Python 3.10, 3.11, and 3.12.
- Automated GitHub Actions workflow passing without warnings.
