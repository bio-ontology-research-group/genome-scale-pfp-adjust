#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=eval_rerun
#SBATCH --array=0-2
#SBATCH --output=%x_%a.log
#
# For one timeset predictor (CC): build the two naive-filter baselines, then run
# a fresh per-organism CAFA + taxon-consistency eval on four settings:
#   predictions (raw) | optimized (CP-SAT, taxon-guided) | naive_prop | naive_noprop
# Feeds: P0 naive-filter comparison (satisfiable %, flip agreement) and
#        P1 paired Wilcoxon/bootstrap (predictions vs optimized).
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):${GAEF_DIR}:${PYTHONPATH:-}"

PREDS=(deepgo-se seq_sim_timeset sprof-go)
THRS_ADJ=(0.3 0.001 0.001)       # tau the predictor was adjusted at (for naive filter)
P=${PREDS[$SLURM_ARRAY_TASK_ID]}
TADJ=${THRS_ADJ[$SLURM_ARRAY_TASK_ID]}
BASE=${DATA_DIR}/swissprot_proteomes_folds/${P}_results/cc
ANN=${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic
NAIVE=${OUT_DIR}/naive/${P}/cc
OUT=${OUT_DIR}/eval/${P}/cc
mkdir -p ${NAIVE} ${OUT}

echo "=== naive filters (tau=${TADJ}) ==="
python analysis/naive_filter.py --predictions_dir ${BASE}/predictions \
    --output_dir ${NAIVE}/prop   --threshold ${TADJ}
python analysis/naive_filter.py --predictions_dir ${BASE}/predictions \
    --output_dir ${NAIVE}/noprop --threshold ${TADJ} --no-propagate

# settings -> source dir (CAFA metrics use eval threshold 0.3 throughout, as in slurm/eval)
declare -A SRC=( [predictions]=${BASE}/predictions [optimized]=${BASE}/optimized [naive_prop]=${NAIVE}/prop [naive_noprop]=${NAIVE}/noprop )
for S in predictions optimized naive_prop naive_noprop; do
    echo "=== eval ${P} ${S} ==="
    python pipeline/evaluate_directory.py \
        --predictions_dir ${SRC[$S]} \
        --annotations_dir ${ANN} \
        --output_dir ${OUT}/${S} \
        --threshold 0.3 --subontology cc \
        --GAEF_dir ${GAEF_DIR} --go_file ${GO_OBO} --num_workers 8 2>&1 | tail -3
done
echo "DONE ${P}"
