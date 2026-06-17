# Reproducibility map: paper claim → producing code

This file maps each table, figure, and headline number in the manuscript to the
script that produces it, per *Bioinformatics* §3. Per-organism intermediate
outputs (predictions, adjusted annotations, per-organism evaluation JSONs) live
under `${DATA_DIR}/swissprot_proteomes_folds/` (see top-level `README.md`);
small aggregated result files for the post-review experiments are committed
under `analysis/revision_results/`.

## Main paper

| Claim / table | Producing code |
|---|---|
| Table 2 — flip rates by predictor × sub-ontology | `slurm/adjust/*.slurm` (`pipeline/run_adjustment_pipeline.py`); rates via `analysis/flip_rate.py` |
| Table 3 — GAEF metrics after adjustment (100% taxon / complex coherence) | `slurm/eval/gaef_timeset_methods.slurm` → `pipeline/evaluate_directory.py` → `pipeline/genome_scores_evaluator.py` |
| §2.4 τ-sensitivity sweep (flip rates across τ) | `slurm/tau_sweep.slurm` (CC), `slurm/tau_sweep_bpmf.slurm` (BP/MF) + `analysis/flip_rate.py` → `revision_results/tau_sweep/` |
| §2.5 naive-filter baseline (satisfiable %) | `analysis/naive_filter.py` + `slurm/eval_rerun.slurm` (`satisfiable` field of `genome_scores_evaluator`) |
| §2.5 paired Wilcoxon / bootstrap (per-organism CAFA) | `slurm/eval_rerun.slurm` / `slurm/eval_all.slurm` → `analysis/paired_stats.py` → `revision_results/wilcoxon/` |
| §2.4 cost-function ablation (margin / uniform / IC) | `analysis/compute_go_ic.py` + `slurm/cost_ablation.slurm` (driven by `PFP_COST_MODE`/`PFP_IC_FILE` in `taxon_consistency/adjust_ortools.py`) → `revision_results/cost_ablation/` |
| §3.1 / Table S8 — GAEF on ground-truth proteomes (6.2%, IC depth) | `slurm/eval/gaef_annotations.slurm` → `genome_scores_evaluator` |
| §3.1 / Table S7 — top-20 curated taxon violations | `analysis/top_violated_pairs.py` → `revision_results/top20_violated_pairs.tsv` |
| §3.4 / Table S9 — per-protein CAFA on modified subset | `pipeline/evaluate_directory.py` (modified-protein subset) |
| §3.6 — Stage-2 biological validity vs Complex Portal | `analysis/complex_validation/` (`download_complex_dbs.sh` → `build_subunit_index.py` → `validate_stage2.py`); broadened run via `slurm/stage2_validation.slurm` → `revision_results/complex_validation/` |
| §3.6 — scalability ablation (5,380× reduction, 12.9 s) | `complex_coherence/heuristic_ablation.py` (`slurm/ablation/heuristic.slurm`) |

## Supplement

| Claim / table | Producing code |
|---|---|
| Tables S1–S2 — per-protein CAFA (timeset / heldout, full) | `pipeline/evaluate_directory.py` |
| Tables S3–S4 — CP-SAT solver runtimes | instrumentation in `taxon_consistency/` and `complex_coherence/adjust_ortools.py` |
| Tables S5–S6 — reduction & heuristic ablation | `complex_coherence/heuristic_ablation.py` |
| §S5 — exactness proofs for H1/H2 | analytic; smoke-checked by `tests/test_solver_smoke.py` |
| §S4 / Fig S1–S2 — joint taxon recovery from contigs | `analysis/metagenomics_*.py` |

## Post-review experiment scripts (added during revision)

- `analysis/flip_rate.py` — Stage-1 flip rate at a given τ.
- `analysis/naive_filter.py` — naive taxon-filter baseline (± GO-hierarchy propagation).
- `analysis/paired_stats.py` — paired Wilcoxon + organism-level bootstrap CIs.
- `analysis/top_violated_pairs.py` — rank curated taxon-constraint violations.
- `analysis/compute_go_ic.py` — GO-term information content for the IC-weighted cost.
- `analysis/complex_validation/validate_stage2.py` — score Stage-2 repairs vs Complex Portal.
- `slurm/{tau_sweep,tau_sweep_bpmf,eval_rerun,eval_all,cost_ablation,stage2_validation,mlp_retrain}.slurm`.
