#!/bin/bash
#SBATCH --job-name=uniprot_predict_deepgo_se
#SBATCH --output=logs/uniprot_predict_deepgo_se_%j.out
#SBATCH --error=logs/uniprot_predict_deepgo_se_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --constraint=[v100|rtx2080ti|gtx1080ti]
#SBATCH --mail-type=ALL

# Usage:
#   sbatch slurm/train/predict_deepgo_se.slurm <TASK>     # TASK ∈ {cc,mf,bp}
# Runs the pretrained DeepGO-SE ensemble on the timeset test organisms.
# Writes ${DATA_DIR}/swissprot_proteomes_folds/deepgo-se_results/<TASK>/predictions/,
# consumed by slurm/adjust/deepgo_se.slurm.

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

# DeepGO-SE pretrained ensemble checkpoints. Override with DEEPGOZERO_MODELS_DIR if needed.
# Default location: ${DEEPGO_SE_REPO}/data/deepgo_se (matches upstream layout).
DEEPGOZERO_MODELS_DIR="${DEEPGOZERO_MODELS_DIR:-${DEEPGO_SE_REPO}/data/deepgo_se}"

python baselines/deepgo_se/predict_deepgo-se.py \
    --subontology ${TASK} \
    --data-root ${DEEPGO_SE_REPO}/data \
    --models-dir ${DEEPGOZERO_MODELS_DIR} \
    --esm2-embeddings-dir ${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings \
    --test-organisms-file splits/timeset/test_organisms.txt \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/deepgo-se_results/${TASK} \
    --go-obo data/go-basic.obo \
    --device cuda:0 \
    --num-workers 8
