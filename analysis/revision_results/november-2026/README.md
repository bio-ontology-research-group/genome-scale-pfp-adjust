# Revised results, November 2026

Code and results of the experiments added after commit `c260eef` for the
revised manuscript (BIOINF-2026-2454). The files under
`analysis/revision_results/` outside this folder are historical.

## paper-tables/

The generated tables and result paragraphs included in the revised manuscript
and its supplement, as TeX and TSV. `generation-manifest.json` records the
inputs and outputs of the table generator for the corrected baselines. Seven
tables were edited or replaced after generation, so their checksums differ from the
manifest: `supp-lineage.tex` (now the lineage objective, see `lineage-data/`),
`supp-accuracy.tex`, `supp-accuracy.tsv`, `supp-comparators.tex`,
`supp-membership.tex`, `supp-solvers.tex` and `supp-threshold-metrics.tex`. The
membership and solver tables received the rows of the evidence-gated Stage 2
(from `reruns/rerun-gated/tables/`); their TSV versions do not contain them.

## lineage-data/

Per-organism lineage-inference results for the repair cost and the lineage
objective; see [`lineage-data/README.md`](lineage-data/README.md).

## reruns/

Scripts and outputs of the follow-up runs. The scripts import the solver and
evaluation code of this repository at `c260eef` (set `PFP_CODE` to a checkout)
and were run on the KAUST IBEX cluster and a workstation; absolute paths in the
scripts and result files refer to those systems. `common.env` lists the
directories used.

| Path | Content |
|---|---|
| `rerun-notaxon/accepted/` | Result bundle of the corrected baselines: per-predictor summaries and the tables generated from them |
| `rerun-notaxon/` | Repair without the organism taxon in two steps: lineage inference with the selected lineage objective, then repair with the inferred taxa fixed (`scripts/`, `summaries/`, `tables/`) |
| `rerun-gated/` | Evidence-gated Stage 2, the default complex-coherence mode of the revised manuscript. `code/revision_audit_run_gated.py` is `analysis/revision_audit_run.py` with the option `--stage2-gated`; `membership/` compares all three modes with curated complex membership |
| `investigate-lineage/` | Comparison of lineage objectives and the selection rule |
| `investigate-complex/` | Stage 2 repair variants and their comparison with Complex Portal and UniProt annotations |
| `investigate-aupr/` | Check of the AUPR computation of the original submission against the corrected evaluator |
| `contig_lineage.py`, `plot_contigs.py`, `results/contigs*/` | Lineage inference on contigs |
| `stage2_ablation.py`, `results/ablation/` | Stage 2 preprocessing ablation |
| `gt_taxon_audit.py`, `results/closure-audit.json` | Closure-based audit of taxon-constraint violations in the reference proteomes |
| `runtime_termcentric.py`, `results/runtime-termcentric.json` | Solver runtime and term-centric analysis of flips |

Job logs and the downloaded FunTaxIS-lite constraint lists
(<https://www.medcomp.medicina.unipd.it/funtaxis-lite/>) are not included.
