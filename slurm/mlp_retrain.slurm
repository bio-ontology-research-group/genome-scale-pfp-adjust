#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=14:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mlp_retrain
#SBATCH --array=0-2
#SBATCH --output=%x_%a.log
#
# Re-train MLP-ESM2 on the heldout split under the fixed RNG seed (--seed 42),
# one sub-ontology per array task, writing to OUT_DIR so Asaad's data is
# untouched. Predictions then feed the adjust+eval reproducibility check.
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):${GAEF_DIR}:${PYTHONPATH:-}"

TASKS=(cc mf bp)
TASK=${TASKS[$SLURM_ARRAY_TASK_ID]}
D=${DATA_DIR}/swissprot_proteomes_folds
OUT=${OUT_DIR}/mlp_seeded/${TASK}
mkdir -p ${OUT}

echo "=== retrain MLP-ESM2 ${TASK} (seed 42) ==="
python baselines/mlp_esm2/train_mlp_heldout.py \
    --subontology ${TASK} \
    --train-file splits/heldout/train_proteins.tsv \
    --val-file   splits/heldout/val_proteins.tsv \
    --test-organisms-file splits/heldout/test_organisms.txt \
    --esm2-embeddings-dir ${D}/esm2_embeddings \
    --annotations-go-basic-dir ${D}/annotations-go-basic \
    --output-dir ${OUT} \
    --device cpu --sparsity-threshold 0.009 --num-workers 8 --seed 42
echo "DONE mlp_retrain ${TASK}"
