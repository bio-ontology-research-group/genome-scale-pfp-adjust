#!/bin/bash
#SBATCH --job-name=uniprot_seq_sim
#SBATCH --output=logs/uniprot_seq_sim_%j.out
#SBATCH --error=logs/uniprot_seq_sim_%j.err
#SBATCH --time=20:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL

# Load modules and activate environment
[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

# Create logs directory
mkdir -p logs

python baselines/seq_sim/predict_seq_similarity.py \
	--subontology all \
	--annotations-dir ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \
	--output-dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results \
	--train-organisms-file splits/heldout/train_organisms.txt \
	--test-organisms-file splits/heldout/test_organisms.txt;
