# Multiplex panel

!!! warning "Computational design phase"
    Panels are **designed**, not experimentally validated. No multiplex
    cross-reactivity has been measured.

## Goal

Combine guide-RNA / primer designs across the priority resistance-gene
families into a **single point-of-care panel** compatible with the
paper-strip assay architecture under design.

## Optimization

`multiplex_panel.py` selects a subset of guides and primers that:

- Cover the prioritized clinical phenotypes (see [library](../data/library.md)).
- Minimize predicted cross-reactivity within the panel.
- Respect amplicon-size and reaction-compatibility constraints for RPA.
- Balance per-gene representation against panel size.

## Clinical interpretation

`clinical_interpreter.py` produces a natural-language summary of a panel
result mapping detected guides to drug class and resistance mechanism via
the [CARD enrichment](scoring.md#card-enrichment) layer.

## Future work

- Cross-reactivity mapping at the wet-lab stage — see
  [Roadmap — Phase 2](../project/roadmap.md#phase-2--panel-expansion).
- Paper-strip / lateral flow prototype.

## Sources

- `multiplex_panel.py` — panel optimization.
- `clinical_interpreter.py` — natural-language interpretation.
