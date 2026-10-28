# genome-scale-pfp-adjust

Post-hoc CP-SAT constraint optimization for repairing genome-scale protein function predictions.

Companion code for the paper *Genome-scale protein function adjustment using constraint optimization*. Implements a two-stage CP-SAT solver (taxon consistency, then complex coherence) that minimally adjusts genome-scale GO-term predictions to satisfy organism-level biological constraints.

## What the solver does

Given per-protein GO predictions from any upstream method, the solver flips the minimum-cost set of annotations so that:

1. **Stage 1 — taxon consistency.** No GO term in the predicted set violates the organism's NCBI Taxonomy lineage (e.g. no *photosynthesis* in an animal proteome). Inferred from `only_in_taxon` / `never_in_taxon` constraints in the Gene Ontology.

2. **Stage 2 — complex coherence.** No obligate heteromeric complex GO term is annotated to exactly one protein in the proteome. Singletons are either demoted from the lone protein or promoted onto a plausible partner.

Both stages are encoded as integer programs and solved with OR-Tools CP-SAT. Stage 2 uses three reductions: incoherent-only filter, sparse variable set, and top-k participation.

## Demo

An executable Showboat demo is available in [`demos/empty_quarter_ibex/README.md`](demos/empty_quarter_ibex/README.md). It fetches a real Empty Quarter genome annotation from the public reproducibility bundle at `https://bio2vec.net/data/genoadjust/rh04_63_rh04_bacillus_spizizenii_genoadjust_demo/` (`rh04` / `63_rh04`, *Bacillus spizizenii*, `NCBITaxon_96241`), converts PGAP GO annotations into the solver's prediction format, and verifies both taxon-consistency and complex-coherence repairs.

## External dependencies

Three baselines and one evaluator load code from sibling repositories. Clone them and point `config.env` at the clones:

- `GAEF_DIR` — `https://github.com/bio-ontology-research-group/GAEF` (only `GAEF.completeness` and `GAEF.utils` are imported from upstream; the divergent modules — `coherence`, `taxon_consistency`, `complex_classifier`, and the `Ontology` IC methods — live under `gaef_patches/` in this repo, and the constraint files under `data/constraints/`).
- `DEEPGO_SE_REPO` — `https://github.com/bio-ontology-research-group/deepgo2`
- `SPROF_GO_REPO` — `https://github.com/biomed-AI/SPROF-GO`

## Quickstart

```
# 1. Environment
conda env create -f environment.yml
conda activate pfp

# 2. Local config
cp config.env.example config.env
$EDITOR config.env                              # set DATA_DIR, CONDA_*, GAEF_DIR, …
source config.env

# 3. Fetch large constraint files
./data_prep/download_data.sh                    # go-basic.obo, ncbitaxon.obo

# 4. Tests (offline, CPU-only, <30 s) — also run in GitHub Actions CI
python data/constraints/verify_checksums.py   # pinned constraint files
python tests/test_data_loaders.py
python tests/test_solver_smoke.py
python tests/test_golden_e2e.py               # 50-protein organism, diffs vs golden output
# (these four need only `pip install -r requirements-test.txt`, not the full conda env)

# 5. Full pipeline
# 5a. Per-protein ESM2 mean embeddings (needed by MLP-ESM2 + DeepGO-SE)
ESM=$(sbatch --parsable slurm/train/extract_esm2.slurm)

# 5b. Upstream predictions (one job per baseline+split)
SS_H=$(sbatch --parsable slurm/train/predict_seq_sim_heldout.slurm)
SS_T=$(sbatch --parsable slurm/train/predict_seq_sim_timeset.slurm)
for T in cc mf bp; do
    sbatch --dependency=afterok:${ESM} slurm/train/train_mlp_heldout.slurm "${T}"
    sbatch --dependency=afterok:${ESM} slurm/train/predict_deepgo_se.slurm "${T}"
done
SPROF=$(sbatch --dependency=afterok:${ESM} --parsable slurm/train/predict_sprof_go.slurm)

# 5c. Adjustment
sbatch --dependency=afterok:${SS_T} slurm/adjust/seq_sim_timeset.slurm
sbatch --dependency=afterok:${SS_H} slurm/adjust/seq_sim_heldout.slurm
sbatch                              slurm/adjust/mlp_heldout.slurm
sbatch                              slurm/adjust/deepgo_se.slurm
for T in cc mf bp; do
    sbatch --dependency=afterok:${SPROF} slurm/adjust/sprof_go.slurm "${T}"
done

# 5d. Evaluation
sbatch slurm/eval/gaef_annotations.slurm
sbatch slurm/eval/gaef_timeset_methods.slurm

# Optional ablation (small/medium/large taxa)
for TX in 272844 759272 10116; do sbatch slurm/ablation/heuristic.slurm "${TX}"; done
```

## Repository layout

```
taxon_consistency/        Stage 1 CP-SAT solver
complex_coherence/        Stage 2 CP-SAT solver + heuristic ablation
pipeline/                 Adjustment orchestrator (chains both stages) + evaluation scripts
gaef_patches/             Local extensions to upstream GAEF (see "External dependencies")
deepgo/                   Utility package (Ontology, MLPBlock, …)
baselines/
  seq_sim/                Sequence-similarity baseline (DIAMOND BLAST)
  mlp_esm2/               MLP trained on ESM2 embeddings.
  deepgo_se/              Wrapper for DeepGO-SE (external repo)
  sprof_go/               Wrapper for SPROF-GO (external repo)
data_prep/                Dataset/split construction, OBO parsers, dataset builder, download_data.sh
data/                     Bundled small constraint TSVs (large ones downloaded)
  constraints/            GAEF constraint files (vendored — see gaef_patches/)
splits/{timeset,heldout}/ Train/test organism splits used in the paper
analysis/                 Supplementary & revision analyses (term-centric, metagenomics, flip ranking, baselines)
reproducibility/          Single-command reproduction of all paper & supplementary tables and figures
slurm/{adjust,eval,train,ablation}/
                          Cluster job scripts, all templated via config.env
tests/                    CPU-only smoke tests + fixtures
```

## `DATA_DIR` layout

`DATA_DIR` is the single working-data root referenced by every SLURM script and
CLI default. You populate `uniprot_reference_proteomes/` yourself (https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/); everything under `swissprot_proteomes_folds/` is produced by
`data_prep/`, the training/prediction jobs, and the adjustment stages.

```
${DATA_DIR}/
├── uniprot_reference_proteomes/        # UniProt reference proteomes as downloaded from the UniProt FTP site
└── swissprot_proteomes_folds/          # all generated artifacts live here
    ├── annotations-go-basic/           # ground-truth GO annotations (data_prep/)
    ├── esm2_embeddings/                # per-protein ESM2 mean embeddings (train/extract_esm2.slurm)
    ├── seq_sim_{timeset,heldout}_results/{cc,mf,bp}/
    │   ├── predictions/                # raw upstream predictions
    │   ├── optimized/                  # after Stage 1 (taxon) + Stage 2 (complex, cc only)
    │   └── optimized_without_taxons/   # after Stage 1 and 2, but without providing taxon information
    ├── mlp_heldout_results/{cc,mf,bp}/<same subdirs>
    ├── deepgo-se_results/{cc,mf,bp}/<same subdirs>
    └── sprof-go_results/{cc,mf,bp}/<same subdirs>
```

SPROF-GO additionally uses `${DATA_DIR}/sprof_go/sprof_go_work/` as an
intermediate working directory (configurable via `WORK_DIR=` on the sbatch
line).

## Reproducing the paper

All tables and figures in the main manuscript (Tables 1–2) and Supplementary Material (Tables S1–S15, Figures S1–S3) can be reproduced from [`reproducibility/`](reproducibility/README.md):

```bash
# Instant offline regeneration of all tables & figures from committed summary artifacts (< 5 s):
uv run python reproducibility/reproduce_all_tables_and_figures.py

# Full end-to-end cluster reproduction (SLURM):
bash reproducibility/run_all_slurm.sh
```

See [`reproducibility/README.md`](reproducibility/README.md) for the complete table-by-table mapping:

| Paper claim | Script |
|---|---|
| Stage 1 & Stage 2 flip rates (Table 1) | `analysis/rev2_tau_and_flips.py` (`slurm/rev2_table1_flips.slurm`) |
| GAEF after adjustment (Table 2; Supp. Table S13) | `analysis/rev2_stage2_and_baselines.py` (`slurm/rev2_stage2_and_baselines.slurm`) |
| Per-protein CAFA metrics (Supp. Tables S1, S2, S9, S13) | `pipeline/evaluate_directory.py` & `analysis/rev2_stage2_and_baselines.py` |
| Metagenomics (Supp. Figures S1–S2) | `analysis/metagenomics_*.py` |
| Heuristic ablation (Supp. Section S6, Tables S5–S6) | `slurm/ablation/heuristic.slurm` & `analysis/ablation_var_counts.py` |
| Threshold sensitivity & cost ablation (Supp. Tables S10–S11) | `analysis/rev2_tau_and_flips.py` (`slurm/rev2_tau_sweep.slurm`) |
| 7-rank taxonomic lineage inference (Supp. Table S12) | `analysis/evaluate_lineage_inference.py` (`slurm/rev2_lineage_and_funtaxis.slurm`) |
| Baselines, FunTaxIS-lite, Soft SBR & Complex Portal (Supp. Tables S13–S15) | `analysis/rev2_stage2_and_baselines.py` & `analysis/funtaxis_lite_filter.py` |

## Citation

> Mohammedsaleh, A. and Hoehndorf, R. (2026). Genome-scale protein function adjustment using constraint optimization. (under preparation).
