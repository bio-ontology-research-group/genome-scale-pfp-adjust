#!/bin/bash
#SBATCH --job-name=uniprot_train_mlp_heldout
#SBATCH --output=logs/uniprot_train_mlp_heldout_%j.out
#SBATCH --error=logs/uniprot_train_mlp_heldout_%j.err
#SBATCH --time=8:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --constraint=[v100|rtx2080ti|gtx1080ti]
#SBATCH --mail-type=ALL

# Usage:
#   sbatch slurm/train/train_mlp_heldout.slurm <TASK>     # TASK ∈ {cc,mf,bp}
# Writes ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/<TASK>/predictions/,
# consumed by slurm/adjust/mlp_heldout.slurm.

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

TASK=$1
if [ -z "${TASK}" ]; then
    echo "Error: TASK is not set"
    exit 1
fi
if [ "${TASK}" != "mf" ] && [ "${TASK}" != "bp" ] && [ "${TASK}" != "cc" ]; then
    echo "Error: TASK must be mf, bp, or cc"
    exit 1
fi

echo "Training MLP-ESM2 on heldout split for ${TASK}..."

python baselines/mlp_esm2/train_mlp_heldout.py \
    --subontology ${TASK} \
    --train-file splits/heldout/train_proteins.tsv \
    --val-file   splits/heldout/val_proteins.tsv \
    --test-organisms-file splits/heldout/test_organisms.txt \
    --esm2-embeddings-dir ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings \
    --annotations-go-basic-dir ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/${TASK} \
    --device cuda:0 \
    --sparsity-threshold 0.009 \
    --num-workers 16
