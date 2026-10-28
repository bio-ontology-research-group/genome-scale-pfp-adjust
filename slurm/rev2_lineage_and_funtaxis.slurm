#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=rev2_lin_ft
#SBATCH --array=0-1
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/rev2/logs/%x_%a.log

set -euo pipefail
cd /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/code
source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):$(pwd)/pipeline:$(dirname ${GAEF_DIR}):${GAEF_DIR}:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

mkdir -p ${OUT_DIR}/rev2/logs ${OUT_DIR}/rev2/funtaxis_lite

TS_FILE="${OUT_DIR}/rev2/timeset_50_taxons.txt"
if [ ! -f "${TS_FILE}" ]; then
    ls ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions/predictions_fold_*_taxon_*.tsv \
        | sed -E 's/.*_taxon_([0-9]+)\.tsv/\1/' | sort -u > "${TS_FILE}.tmp.$$"
    mv "${TS_FILE}.tmp.$$" "${TS_FILE}"
fi

if [ "${SLURM_ARRAY_TASK_ID}" = "0" ]; then
    echo "=== Running evaluate_lineage_inference.py ==="
    python analysis/evaluate_lineage_inference.py \
        --data_dir "${DATA_DIR}" \
        --code_dir "$(pwd)" \
        --timeset_taxons_file "${TS_FILE}" \
        --output_json "${OUT_DIR}/rev2/lineage_inference.json" \
        --output_tsv  "${OUT_DIR}/rev2/lineage_inference.tsv"
else
    if [ ! -s "${OUT_DIR}/rev2/funtaxis_lite_summary.json" ]; then
        echo "=== Running funtaxis_lite_filter.py ==="
        python analysis/funtaxis_lite_filter.py \
            --data_dir "${DATA_DIR}" \
            --code_dir "$(pwd)" \
            --timeset_taxons_file "${TS_FILE}" \
            --out_root "${OUT_DIR}/rev2/funtaxis_lite" \
            --summary_json "${OUT_DIR}/rev2/funtaxis_lite_summary.json" \
            --cutoff 500
    fi
fi

echo "DONE rev2_lin_ft task ${SLURM_ARRAY_TASK_ID}"
