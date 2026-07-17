# MOSAIC Layer 0 — Biological Entities

Cells

↓

Biological Entities

## Purpose

Layer 0 defines the biological entities that enter the MOSAIC framework before any graph construction or network analysis takes place.

MOSAIC is intentionally modality-agnostic. Any technology capable of identifying biological objects with spatial coordinates can serve as input, provided that the required metadata are available.

## Core Entity: Cell

Each detected cell represents the fundamental unit of the graph.

Required attributes include:

* Cell ID
* X coordinate
* Y coordinate
* Cell phenotype (or cell type)
* Sample or slide identifier

## Optional Attributes

Additional metadata can be incorporated when available:

* Tissue compartment
* Disease type
* Patient identifier
* Clinical outcome
* Imaging platform
* Marker intensities
* Gene expression
* Protein abundance
* Mutation status
* Cell quality metrics
* Time point
* Experimental batch

## Design Principle

Layer 0 does not perform biological interpretation.

Its purpose is to establish a standardized representation of biological entities that can be transformed into spatial graphs in subsequent layers.

## Output

A harmonized collection of biologically annotated spatial objects that serve as the input for Layer 1: Spatial Graph Construction.

