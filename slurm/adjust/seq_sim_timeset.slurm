#!/bin/bash
#SBATCH --job-name=uniprot_timeset_seq_sim_adjust
#SBATCH --output=logs/uniprot_timeset_seq_sim_adjust_%j.out
#SBATCH --error=logs/uniprot_timeset_seq_sim_adjust_%j.err
#SBATCH --time=6:00:00
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

BASE=${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results

# Stage 1, taxon-guided (writes optimized/)
for TASK in cc mf bp; do
    echo "Seq-Sim ${TASK} | Stage 1 with taxon ID"
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE}/${TASK}/predictions \
        --optimized_dir   ${BASE}/${TASK}/optimized \
        --optimized --default_threshold 0.001 --provide_taxon_id
done

# Stage 1, taxon-unknown (writes optimized_without_taxons/)
for TASK in cc mf bp; do
    echo "Seq-Sim ${TASK} | Stage 1 without taxon ID"
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE}/${TASK}/predictions \
        --optimized_dir   ${BASE}/${TASK}/optimized_without_taxons \
        --optimized --default_threshold 0.001
done

# Stage 2 (complex coherence) on top of Stage 1 taxon-unknown CC only
echo "Seq-Sim CC | Stage 1+2 (complex coherence, top_k=50, without taxon ID)"
python pipeline/run_adjustment_pipeline.py \
    --predictions_dir ${BASE}/cc/predictions \
    --optimized_dir   ${BASE}/cc/optimized_without_taxons_complex_coherence \
    --optimized --complex_coherence --top_k 50 --default_threshold 0.001

