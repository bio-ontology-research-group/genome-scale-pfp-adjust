#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mlp_reproduce
#SBATCH --output=%x_%j.log
#
# Faithful, seeded reproduction of the MLP-ESM2 baseline (paper Section 2.2).
# 1. Regenerate the per-protein train/val split deterministically from the
#    preserved organism lists (test set held fixed; MMseqs2 single-threaded;
#    cluster split seeded at 42).
# 2. Train MLP-ESM2 for cc/mf/bp under --seed 42 and predict on the held-out
#    test organisms.
# Outputs go to OUT_DIR so Asaad's read-only data is untouched. Adjust+eval is a
# separate step (slurm/mlp_reproduce_eval.slurm) once predictions exist.
set -euo pipefail
[ -f config.env ] && source config.env
# Load the mmseqs module FIRST so that conda activate (which prepends the env bin)
# wins on PATH; otherwise the module's own python shadows the conda python.
module load mmseqs2/14.7 2>/dev/null || true
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
PY=${CONDA_BASE}/envs/${CONDA_ENV}/bin/python
export PYTHONPATH="$(pwd):${GAEF_DIR}:${PYTHONPATH:-}"
echo "python: $(command -v ${PY}) ; mmseqs: $(command -v mmseqs || echo MISSING)"
${PY} -c "import pandas; print('pandas', pandas.__version__)"

SPLIT=splits/heldout
D=${DATA_DIR}/swissprot_proteomes_folds
OUT=${OUT_DIR}/mlp_seeded
mkdir -p ${OUT}

echo "=== Step 1: regenerate train/val protein split (deterministic, seed 42) ==="
if [ ! -s ${SPLIT}/train_proteins.tsv ] || [ ! -s ${SPLIT}/val_proteins.tsv ]; then
    ${PY} data_prep/regen_heldout_proteins.py \
        --split_dir ${SPLIT} \
        --proteomes_ids_file data/uniprot_proteomes_ids.tsv \
        --annotations_dir ${D}/annotations-go-basic \
        --main_proteomes_dir ${DATA_DIR}/uniprot_reference_proteomes \
        --threads 1 --seed 42
else
    echo "train/val TSVs already present, skipping regen"
fi

# The MLP output vocabulary (GO terms per sub-ontology) lives in the deepgo2-data
# bundle; the training script expects it at data/<ont>/terms.pkl. Symlink it.
TERMS_SRC=${DATA_DIR}/deepgo2-data/data
for TASK in cc mf bp; do
    mkdir -p data/${TASK}
    [ -e data/${TASK}/terms.pkl ] || ln -sf ${TERMS_SRC}/${TASK}/terms.pkl data/${TASK}/terms.pkl
done

for TASK in cc mf bp; do
    echo "=== Step 2: train MLP-ESM2 ${TASK} (seed 42) ==="
    mkdir -p ${OUT}/${TASK}
    ${PY} baselines/mlp_esm2/train_mlp_heldout.py \
        --subontology ${TASK} \
        --train-file ${SPLIT}/train_proteins.tsv \
        --val-file   ${SPLIT}/val_proteins.tsv \
        --test-organisms-file ${SPLIT}/test_organisms.txt \
        --esm2-embeddings-dir ${D}/esm2_embeddings \
        --annotations-go-basic-dir ${D}/annotations-go-basic \
        --output-dir ${OUT}/${TASK} \
        --device cpu --sparsity-threshold 0.009 --num-workers 16 --seed 42
done
echo "DONE mlp_reproduce: predictions under ${OUT}/{cc,mf,bp}"
