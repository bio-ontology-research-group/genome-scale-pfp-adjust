#!/bin/bash
#SBATCH --job-name=eval_annots_gaef
#SBATCH --output=logs/eval_annots_gaef_%j.out
#SBATCH --error=logs/eval_annots_gaef_%j.err
#SBATCH --time=12:00:00
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

BASE_DIR=${DATA_DIR}/swissprot_proteomes_folds

python pipeline/evaluate_annotations_gaef.py \
	--annotations_dir ${BASE_DIR}/annotations-go-basic \
	--output_dir ${BASE_DIR}/annotations-go-basic/gaef_evaluations \
	--go_file data/go-basic.obo \
	--GAEF_dir ${GAEF_DIR} \
	--num_workers 16;
