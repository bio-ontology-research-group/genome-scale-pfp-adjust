#!/bin/bash
#SBATCH --job-name=uniprot_extract_esm2
#SBATCH --output=logs/uniprot_extract_esm2_%j.out
#SBATCH --error=logs/uniprot_extract_esm2_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --constraint=[v100|rtx2080ti|gtx1080ti]
#SBATCH --mail-type=ALL

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

# Extract ESM2 mean embeddings for every proteome in data/uniprot_proteomes_ids.tsv.
# Outputs ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings/<taxon_id>/esm2.pkl,
# consumed by train_mlp_heldout.slurm and predict_deepgo_se.slurm.

python baselines/mlp_esm2/extract_esm2_embeddings.py \
    --proteomes-ids-file data/uniprot_proteomes_ids.tsv \
    --input-dir ${DATA_DIR}/uniprot_reference_proteomes \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings \
    --device cuda:0 \
    --batch-size 4096
