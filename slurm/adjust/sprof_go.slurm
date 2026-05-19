#!/bin/bash
#SBATCH --job-name=uniprot_sprof-go_adjust
#SBATCH --output=logs/uniprot_sprof-go_adjust_%j.out
#SBATCH --error=logs/uniprot_sprof-go_adjust_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL

# Load modules and activate environment
[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Create logs directory
mkdir -p logs

# get TASK from command line
TASK=$1

if [ -z "${TASK}" ]; then
    echo "Error: TASK is not set"
    exit 1
fi

if [ "${TASK}" != "cc" ] && [ "${TASK}" != "mf" ] && [ "${TASK}" != "bp" ]; then
    echo "Error: TASK must be cc, mf, or bp"
    exit 1
fi

BASE=${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results

echo "SPROF-GO ${TASK} | Stage 1 with taxon ID"
python pipeline/run_adjustment_pipeline.py \
    --predictions_dir ${BASE}/${TASK}/predictions \
    --optimized_dir   ${BASE}/${TASK}/optimized \
    --optimized --default_threshold 0.001 --provide_taxon_id

echo "SPROF-GO ${TASK} | Stage 1 without taxon ID"
python pipeline/run_adjustment_pipeline.py \
    --predictions_dir ${BASE}/${TASK}/predictions \
    --optimized_dir   ${BASE}/${TASK}/optimized_without_taxons \
    --optimized --default_threshold 0.001

# Stage 1+2 (complex coherence) on CC only
if [ "${TASK}" = "cc" ]; then
    echo "SPROF-GO CC | Stage 1+2 (complex coherence, top_k=50, without taxon ID)"
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE}/cc/predictions \
        --optimized_dir   ${BASE}/cc/optimized_without_taxons_complex_coherence \
        --optimized --complex_coherence --top_k 50 --default_threshold 0.001
fi

