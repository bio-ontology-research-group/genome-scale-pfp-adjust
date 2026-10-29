# Repository guidance

## Audited revision workflow

Use `reproducibility/run_all_slurm.sh` with an explicit configuration based on `reproducibility/revision.env.example`. The current workflow reruns adjustment and evaluation on fixed predictor inputs; upstream training is separate. Read `reproducibility/README.md` for the complete input and dependency contract.

Historical numbers in `analysis/revision_results/` and `reproducibility/output/` were invalidated by the October 2026 audit. Do not copy them into manuscript tables. Table generation requires the complete accepted result bundle and cross-checks hashes and counts. Source and inputs must remain immutable during a run.

## Architecture

- `gaef_patches/ontology.py` explicitly adapts the ontology. It resolves alternative IDs to primary IDs, drops unknown/obsolete IDs, merges duplicate primary-term scores by maximum, and supports maximum-score ancestry propagation. There is no import-time IC monkeypatch. Inheritance uses `is_a` and `part_of`, not `regulates`.
- Stage 1 (`taxon_consistency/adjust_ortools.py`) repairs genome-level taxon incompatibility. It accepts a known lineage or jointly infers compatible taxon variables. Union membership is disjunctive; union members are not inherently disjoint. Export removes demoted terms and their descendants across every protein and score level.
- Stage 2 (`complex_coherence/adjust_ortools.py`) enforces zero or at least two proteins per operational complex term. `protein_complexes.tsv` supplies homomeric/ambiguous exclusions, not a positive list of obligate heteromers. The vocabulary is derived from `is_a` descendants of GO:0032991. Cardinality satisfaction does not establish subunit identity or interaction.
- Sparse Stage 2 closes its domain over hierarchy and affected complex constraints. Top-k promotion restriction remains a heuristic. Feasible bounded solutions are not proven optima; retain objective bounds and actual search budgets. No solution must raise an error.
- The revised runner (`analysis/revision_audit_run.py`) canonicalizes and ancestor-propagates all inputs before comparisons. Lower-level TSV APIs expect aligned inputs. The strict structural activation rule is `score > tau`.
- Accuracy uses explicit test proteins with non-root reference labels; structural checks use all predicted proteins. Keep the known-taxon check separate from GAEF's pairwise screen. Fixed-corpus IC and continuous-score PR conventions are documented in the supplement and evaluator tests.
- The public orchestrator is `pipeline/run_adjustment_pipeline.py`. Unverified `--skip_existing` reuse is disabled. `--provide_taxon_id` selects known-taxon repair; `--demotion-only` prevents Stage 2 promotions.

## Environment and checks

`GAEF_DIR` identifies the archived upstream GAEF checkout; local divergent modules live in `gaef_patches/`, with reference resources in `data/constraints/`. Upstream predictor wrappers additionally require their external repositories. Do not replace archived resources with current downloads when claiming an identical rerun.

Install `requirements-test.txt`, then run the nine checks listed in `reproducibility/README.md`. They cover loaders/checksums, solvers, the golden fixture, reviewer regressions, metrics, lineage, real shell dependencies and complete-bundle generation. Passing small tests does not replace the real-data run or final numerical reconciliation.

Use project-local scratch or `$TMPDIR`, never `/tmp`.

## Input format

Prediction rows are `protein_id<TAB>GO:0001234|0.87<TAB>...`. Raw filenames are `predictions_fold_<FOLD>_taxon_<TAXON>.tsv`; repaired filenames are `optimized_fold_<FOLD>_taxon_<TAXON>.tsv`. Preserve exact protein IDs and floating-point scores. See each CLI's `--help` for required paths; do not infer the revision's configuration from older exploratory SLURM wrappers.
