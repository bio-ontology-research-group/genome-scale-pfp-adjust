# Lineage inference: per-organism data

Per-organism data behind Supplementary Tables S13 to S15 of the revised
manuscript (joint lineage inference when the organism is not given).

| File | One row per | Content |
|---|---|---|
| `lineage-ranks-repair-cost.tsv` | predictor, organism, sub-ontology and rank (16,464 rows) | true taxon, taxa required by the retained GO terms and solver-assigned taxa, with names; repair cost as objective |
| `lineage-outcomes-repair-cost.tsv` | predictor, organism and sub-ontology (2,352 rows) | deepest correct rank, first divergence, prefix status, compatibility across sub-ontologies; repair cost as objective |
| `lineage-ranks-all-optima-repair-cost.tsv` | predictor, organism, sub-ontology and rank (16,464 rows) | taxa required by the selected repair and by every optimal repair |
| `lineage-outcomes-by-objective.tsv` | objective, predictor, organism and sub-ontology (4,704 rows) | outcome (correct prefix, wrong, unresolved), deepest correct rank, first divergence, and for the lineage objective the inferred taxa per rank |

The sub-ontology `joint` is the solve over CC, MF and BP together. In all
files a taxon counts as inferred only if every optimal solution requires it.

`make_objective_table.py` writes `lineage-outcomes-by-objective.tsv` and panel
(b) of Supplementary Table S15 from `lineage-outcomes-repair-cost.tsv` and the
output of `../reruns/investigate-lineage/lineage_variants.py` (variant
`siglogit`).

`make_lineage_objective_table.py` writes Supplementary Table S13
(`../paper-tables/supp-lineage.tex`), which reports the lineage objective, from
the same records. `supp-lineage-repair-cost.tex` is the earlier version of that
table for the repair cost.
