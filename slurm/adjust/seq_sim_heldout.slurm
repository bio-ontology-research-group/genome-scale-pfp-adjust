#!/bin/bash
#SBATCH --job-name=uniprot_heldout_seq-sim_adjust
#SBATCH --output=logs/uniprot_heldout_seq-sim_adjust_%j.out
#SBATCH --error=logs/uniprot_heldout_seq-sim_adjust_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL

# Load modules and activate environment
[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Create logs directory
mkdir -p logs

for TASK in cc mf bp; do
    echo "Adjusting held-out predictions with sequence similarity for ${TASK}..."
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results/${TASK}/predictions \
        --optimized_dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results/${TASK}/optimized \
        --optimized --default_threshold 0.001 --provide_taxon_id
done

for TASK in cc mf bp; do
    echo "Adjusting held-out predictions with sequence similarity without providing taxon_id for ${TASK}..."
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results/${TASK}/predictions \
        --optimized_dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results/${TASK}/optimized_without_taxons \
        --optimized --default_threshold 0.001
done

