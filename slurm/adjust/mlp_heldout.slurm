#!/bin/bash
#SBATCH --job-name=uniprot_heldout_mlp_adjust
#SBATCH --output=logs/uniprot_heldout_mlp_adjust_%j.out
#SBATCH --error=logs/uniprot_heldout_mlp_adjust_%j.err
#SBATCH --time=8:00:00
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

# Usage: sbatch uniprot_heldout_mlp_adjust.slurm [with-taxons|without-taxons|both]
# Default: both
MODE=${1:-both}

run_without_taxons() {
    for TASK in cc mf bp; do
        echo "Adjusting heldout predictions with MLP without providing taxon_id for ${TASK}..."
        python pipeline/run_adjustment_pipeline.py \
            --predictions_dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/${TASK}/predictions \
            --optimized_dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/${TASK}/optimized_without_taxons \
            --optimized --default_threshold 0.001
    done
}

run_with_taxons() {
    for TASK in cc mf bp; do
        echo "Adjusting heldout predictions with MLP providing taxon_id for ${TASK}..."
        python pipeline/run_adjustment_pipeline.py \
            --predictions_dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/${TASK}/predictions \
            --optimized_dir ${DATA_DIR}/swissprot_proteomes_folds/mlp_heldout_results/${TASK}/optimized \
            --optimized --default_threshold 0.001 --provide_taxon_id
    done
}

case "$MODE" in
    with-taxons)
        run_with_taxons
        ;;
    without-taxons)
        run_without_taxons
        ;;
    both)
        run_without_taxons
        run_with_taxons
        ;;
    *)
        echo "Unknown mode: $MODE. Use: with-taxons, without-taxons, or both"
        exit 1
        ;;
esac
