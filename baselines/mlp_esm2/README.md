# MLP-ESM2 baseline

A 3-layer MLP trained from scratch on mean-pooled ESM2 (`esm2_t36_3B_UR50D`)
embeddings. Trained per-ontology (`cc` / `mf` / `bp`) on the heldout split and
predicted on its 219 test organisms. This is the only baseline trained inside
this repo; the other three use pretrained external models.

## Files

- `extract_esm2_embeddings.py` — one-time GPU pass over every reference
  proteome listed in `data/uniprot_proteomes_ids.tsv`, writing
  `<output-dir>/<taxon_id>/esm2.pkl` (one 2560-d mean embedding per protein).
- `train_mlp_heldout.py` — trains the per-ontology MLP on the heldout
  train/val partition and writes `predictions_fold_01_taxon_<id>.tsv` for
  every test organism, in the format consumed by
  `pipeline/run_adjustment_pipeline.py`.

## Setup

No external repo. Just the conda env from the top-level `environment.yml`
(`fair-esm`, `torch`, `pytorch-cuda=12.1`) and a CUDA-capable GPU.

## Run via SLURM

```
# 1. Per-protein ESM2 embeddings (run once, ~hours)
sbatch slurm/train/extract_esm2.slurm

# 2. Train + predict per ontology (depends on step 1)
for T in cc mf bp; do
    sbatch --dependency=afterok:${ESM} slurm/train/train_mlp_heldout.slurm "${T}"
done
```

The slurm wrappers pass through the required paths from `config.env`
(`DATA_DIR`, `CONDA_ENV`). Output goes to
`${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/<TASK>/`.

## Direct CLI

```
python baselines/mlp_esm2/extract_esm2_embeddings.py \
    --proteomes-ids-file data/uniprot_proteomes_ids.tsv \
    --input-dir  ${DATA_DIR}/uniprot_reference_proteomes \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings

python baselines/mlp_esm2/train_mlp_heldout.py \
    --subontology cc \
    --train-file splits/heldout/train_proteins.tsv \
    --val-file   splits/heldout/val_proteins.tsv \
    --test-organisms-file splits/heldout/test_organisms.txt \
    --esm2-embeddings-dir ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/cc
```

## Reference

Lin, Z. et al. (2023). Evolutionary-scale prediction of atomic-level protein
structure. *Science* 379:1123-1130. (ESM2 backbone)
