#!/bin/bash
#SBATCH --job-name=uniprot_seq_sim_timeset
#SBATCH --output=logs/uniprot_seq_sim_timeset_%j.out
#SBATCH --error=logs/uniprot_seq_sim_timeset_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

# Diamond-BLAST Seq-Sim baseline on the timeset split.
# Writes ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/{cc,mf,bp}/predictions/,
# consumed by slurm/adjust/seq_sim_timeset.slurm.

python baselines/seq_sim/predict_seq_similarity.py \
    --subontology all \
    --annotations-dir ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \
    --output-dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results \
    --train-organisms-file splits/timeset/train_organisms.txt \
    --test-organisms-file splits/timeset/test_organisms.txt
