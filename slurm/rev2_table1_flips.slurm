#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=rev2_flips
#SBATCH --array=0-14
#SBATCH --output=slurm-logs/%x_%A_%a.log

set -euo pipefail
printf "This historical workflow is retired; use reproducibility/run_all_slurm.sh with revision.env.\n" >&2
exit 2
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):$(pwd)/pipeline:$(dirname ${GAEF_DIR}):${GAEF_DIR}:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

mkdir -p ${OUT_DIR}/rev2/logs ${OUT_DIR}/rev2/stage1

TS_FILE="${OUT_DIR}/rev2/timeset_50_taxons.txt"
if [ ! -f "${TS_FILE}" ]; then
    ls ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions/predictions_fold_*_taxon_*.tsv \
        | sed -E 's/.*_taxon_([0-9]+)\.tsv/\1/' | sort -u > "${TS_FILE}.tmp.$$"
    mv "${TS_FILE}.tmp.$$" "${TS_FILE}"
fi

COMBOS=(
    "seq_sim_timeset:cc" "seq_sim_timeset:mf" "seq_sim_timeset:bp"
    "deepgo-se:cc"       "deepgo-se:mf"       "deepgo-se:bp"
    "sprof-go:cc"        "sprof-go:mf"        "sprof-go:bp"
    "seq_sim_heldout:cc" "seq_sim_heldout:mf" "seq_sim_heldout:bp"
    "mlp_heldout:cc"     "mlp_heldout:mf"     "mlp_heldout:bp"
)

SPEC="${COMBOS[$SLURM_ARRAY_TASK_ID]}"
P="${SPEC%%:*}"
O="${SPEC##*:}"
TEMPORAL_ARGS=()
if [[ "$P" != *_heldout ]]; then
    TEMPORAL_ARGS=(--temporal-proteins "$PWD/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv")
fi

if [ ! -s "${OUT_DIR}/rev2/flips_${P}_${O}.json" ]; then
    echo "=== Running table1_flips for ${P} ${O} ==="
    python analysis/rev2_tau_and_flips.py "${TEMPORAL_ARGS[@]}" \
        --data_dir "${DATA_DIR}" \
        --code_dir "$(pwd)" \
        --timeset_taxons_file "${TS_FILE}" \
        --task table1_flips \
        --predictor "${P}" \
        --subontology "${O}" \
        --save_adjusted_root "${OUT_DIR}/rev2/stage1" \
        --output_json "${OUT_DIR}/rev2/flips_${P}_${O}.json"
fi

if [ "${P}" = "deepgo-se" ] && [ "${O}" = "cc" ] && [ ! -s "${OUT_DIR}/rev2/cost_ablation_deepgo-se_cc.json" ]; then
    echo "=== Running cost_ablation for deepgo-se cc ==="
    python analysis/rev2_tau_and_flips.py "${TEMPORAL_ARGS[@]}" \
        --data_dir "${DATA_DIR}" \
        --code_dir "$(pwd)" \
        --timeset_taxons_file "${TS_FILE}" \
        --task cost_ablation \
        --predictor deepgo-se \
        --subontology cc \
        --output_json "${OUT_DIR}/rev2/cost_ablation_deepgo-se_cc.json"
fi

echo "DONE rev2_flips ${P} ${O}"
