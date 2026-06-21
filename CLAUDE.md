# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment and config

All entry points (Python CLIs, SLURM scripts, `data_prep/download_data.sh`) assume `config.env` has been sourced. Copy `config.env.example` → `config.env`, fill in paths, then `source config.env` before running anything. The variables that downstream code reads:

- `DATA_DIR` — root of all working data (predictions, embeddings, intermediate, adjusted outputs). The layout under `${DATA_DIR}/swissprot_proteomes_folds/` is documented in `README.md`.
- `CONDA_BASE`, `CONDA_ENV` — SLURM scripts do `source ${CONDA_BASE}/etc/profile.d/conda.sh && conda activate ${CONDA_ENV}`.
- `GAEF_DIR`, `DEEPGO_SE_REPO`, `SPROF_GO_REPO` — paths to **sibling clones** of three external repos. Wrappers in `baselines/{deepgo_se,sprof_go}/` and `pipeline/evaluate_*` import from these clones; they are not pip-installed.
- `GO_OBO` — defaults to `data/go-basic.obo` (fetched by `data_prep/download_data.sh`).

Conda env: `conda env create -f environment.yml` → `conda activate pfp` (Python 3.10, PyTorch + CUDA 12.1, OR-Tools, DIAMOND, MMseqs2, fair-esm).

## Smoke tests (CPU, offline)

```bash
python tests/test_data_loaders.py
python tests/test_solver_smoke.py
```

`test_solver_smoke.py` exercises both stages end-to-end on tiny fixtures in `tests/fixtures/` and should finish in seconds. There is no pytest harness — tests are plain scripts with `if __name__ == "__main__"` blocks.

## Architecture: two-stage CP-SAT adjustment

The core contribution is a two-stage min-flip integer program over per-protein GO predictions. Understanding the call graph requires reading several files together:

1. **Stage 1 — `taxon_consistency/adjust_ortools.py`.** Operates on one organism at a time. Loads `data/constraints/go_taxon_constraints_updated.tsv` (`only_in_taxon` / `never_in_taxon`), the GO `is_a`/`part_of` hierarchy, and the NCBI Taxonomy hierarchy (`is_a` + `disjoint_from` + `union_of`). Flips the minimum-cost set of GO annotations so no prediction survives that contradicts the organism's lineage. The slow loaders are factored out (`load_constraints`, `load_go_hierarchy`, `load_taxon_hierarchy`) so the orchestrator can call `adjust_per_taxon(...)` in a loop without re-parsing.

2. **Stage 2 — `complex_coherence/adjust_ortools.py`.** Runs only on the CC (cellular component) namespace. Reads `complex_coherence/protein_complexes.tsv` (obligate heteromeric complexes). For each genome, ensures no such GO term is annotated to exactly one protein — either demoting the singleton or promoting a partner. Stage 2 has three knobs that materially change problem size: incoherent-only filter, sparse variable set, and `top_k` participation; pass them via `optimized=True` and `top_k=<N>` to `complex_coherence_adjust`.

3. **Orchestrator — `pipeline/run_adjustment_pipeline.py`.** Globs `predictions_fold_*_taxon_*.tsv` from a predictions directory, sorts taxon/fold pairs for deterministic parallelism, then loops calling Stage 1 (and Stage 2, if `--complex_coherence`) per taxon. Supports `--start_index`/`--end_index` for SLURM array sharding and `--skip_existing` for resumability. `--provide_taxon_id` toggles whether the organism's NCBI taxon is supplied to Stage 1 (the "with-taxon" vs "without-taxon" experiments).

Heuristic baselines for the ablation section live in `complex_coherence/heuristic_ablation.py` (H1/H2/H3), not in the CP-SAT module.

## `gaef_patches/` — runtime monkey-patches of upstream GAEF

`gaef_patches/__init__.py` imports `GAEF.utils.Ontology` from the upstream `GAEF_DIR` clone and **attaches `calculate_ic` / `get_ic` / `get_norm_ic` methods to that class at import time** (idempotently). Anywhere in this repo that does `from gaef_patches import ...`, the upstream `Ontology` class gains IC methods as a side effect. The other files in `gaef_patches/` (`coherence.py`, `taxon_consistency.py`, `complex_classifier.py`) are local divergent ports — import these instead of the upstream equivalents, even though `GAEF.completeness` and `GAEF.utils` are still imported from upstream. If you find code that imports `from GAEF.coherence` or `from GAEF.taxon_consistency`, that is a bug.

The vendored constraint files referenced by these modules live in `data/constraints/`.

## Prediction file format

Every predictor and every stage of the pipeline reads/writes the same TSV format:

```
protein_id<TAB>GO:0001234|0.87<TAB>GO:0005678|0.42<TAB>...
```

Per-organism prediction files are named `predictions_fold_<FOLD>_taxon_<NCBI_TAXON_ID>.tsv`; per-organism adjusted outputs are `optimized_fold_<FOLD>_taxon_<NCBI_TAXON_ID>.tsv`. The orchestrator parses fold + taxon directly out of these filenames — preserve the naming convention.

## SLURM phases

The pipeline is split into four phases (see `slurm/README.md` for the full table):

- `slurm/train/extract_esm2.slurm` — per-protein ESM2 mean embeddings (input to MLP-ESM2 and DeepGO-SE).
- `slurm/train/predict_*.slurm` and `slurm/train/train_mlp_heldout.slurm` — produce raw upstream predictions per baseline + split. The MLP, DeepGO-SE, and SPROF-GO jobs take a `TASK` arg (`cc`/`mf`/`bp`).
- `slurm/adjust/*.slurm` — call `pipeline/run_adjustment_pipeline.py` to run Stage 1 (with and without taxon info) and Stage 2 (CC only).
- `slurm/eval/*.slurm` — GAEF metrics over ground truth and adjusted predictions.

All SLURM scripts `source config.env` at the top; per-script overrides are documented in each script's header.

## Single-test / single-job execution

- One Stage-1 run on a single taxon: invoke `taxon_consistency/adjust_ortools.py` directly with `--predictions_file`, `--constraints_file`, `--go_hierarchy_file`, `--taxon_hierarchy_file`, `--ncbitaxon_hierarchy_file`, `--output_file`, `--threshold` (use `python <file> --help` for the full list). The same pattern applies to `complex_coherence/adjust_ortools.py`.
- One taxon through the orchestrator: pass `--taxon_ids_file <file with one taxon per line>` to `pipeline/run_adjustment_pipeline.py`.
- SLURM array shard locally: pass `--start_index N --end_index M`.
