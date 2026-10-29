#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=rev2_s2_base
#SBATCH --array=0-2
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

PREDS=(seq_sim_timeset deepgo-se sprof-go)
P="${PREDS[$SLURM_ARRAY_TASK_ID]}"
TS_FILE="${OUT_DIR}/rev2/timeset_50_taxons.txt"

echo "=== Running rev2_stage2_and_baselines.py for ${P} ==="
python analysis/rev2_stage2_and_baselines.py \
    --temporal-proteins "$PWD/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv" \
    --data_dir "${DATA_DIR}" \
    --code_dir "$(pwd)" \
    --gaef_dir "${GAEF_DIR}" \
    --timeset_taxons_file "${TS_FILE}" \
    --stage1_root "${OUT_DIR}/rev2/stage1" \
    --funtaxis_root "${OUT_DIR}/rev2/funtaxis_lite" \
    --naive_root "${OUT_DIR}/naive" \
    --predictor "${P}" \
    --output_json "${OUT_DIR}/rev2/stage2_and_baselines_${P}.json"


# Run Complex Portal validation (validate_stage2.py) on both standard and demotion_only first
IDX="${CP_SUBUNIT_INDEX:?Set CP_SUBUNIT_INDEX in config.env}"
TAU="0.001"
if [ "${P}" = "deepgo-se" ]; then
    TAU="0.3"
fi

for MODE in standard demotion_only; do
    for S1L in notax withtax; do
        S1D="${OUT_DIR}/rev2/stage1/${P}/cc/optimized"
        if [ "${S1L}" = "notax" ]; then
            S1D="${OUT_DIR}/rev2/stage1/${P}/cc/optimized_without_taxons"
        fi
        S2D="${OUT_DIR}/rev2/stage1/${P}/cc/stage2_${MODE}_${S1L}"
        if [ -d "${S2D}" ] && [ ! -f "${OUT_DIR}/rev2/cp_val_${P}_${MODE}_${S1L}.json" ]; then
            python analysis/complex_validation/validate_stage2.py \
                --stage1_dir "${S1D}" \
                --stage2_dir "${S2D}" \
                --subunit_index "${IDX}" \
                --protein-id-map "${PROTEIN_ID_MAP:?Set PROTEIN_ID_MAP in config.env}" \
                --threshold "${TAU}" \
                --label "${P}_${MODE}_${S1L}" \
                --output_file "${OUT_DIR}/rev2/cp_val_${P}_${MODE}_${S1L}.json"
        fi
    done
done

# For deepgo-se and sprof-go, also run demotion-only on the 6 Complex Portal reference organisms
if [ "${P}" = "deepgo-se" ] || [ "${P}" = "sprof-go" ]; then
    CP_S1="${DATA_DIR}/swissprot_proteomes_folds/${P}_results/cc/optimized_without_taxons"
    CP_S2_STD="${OUT_DIR}/stage2_validation/${P}/cc/stage2"
    CP_S2_DEM="${OUT_DIR}/rev2/cp6_${P}_demotion_only"
    if [ -d "${CP_S2_STD}" ] && [ ! -f "${OUT_DIR}/rev2/cp6_val_${P}_demotion_only.json" ]; then
        mkdir -p "${CP_S2_DEM}"
        for f in ${CP_S2_STD}/optimized_fold_*_taxon_*.tsv; do
            bn=$(basename "$f")
            if [ -f "${CP_S1}/${bn}" ]; then
                python complex_coherence/adjust_ortools.py \
                    --predictions "${CP_S1}/${bn}" \
                    --complexes complex_coherence/protein_complexes.tsv \
                    --go_hierarchy data/go_hierarchy_cc.tsv \
                    --output "${CP_S2_DEM}/${bn}" \
                    --threshold "${TAU}" \
                    --optimized --top_k 50 --demotion_only
            fi
        done
        python analysis/complex_validation/validate_stage2.py \
            --stage1_dir "${CP_S1}" \
            --stage2_dir "${CP_S2_DEM}" \
            --subunit_index "${IDX}" \
                --protein-id-map "${PROTEIN_ID_MAP:?Set PROTEIN_ID_MAP in config.env}" \
            --threshold "${TAU}" \
            --label "${P}_cp6_demotion_only" \
            --output_file "${OUT_DIR}/rev2/cp6_val_${P}_demotion_only.json"
    fi
fi


echo "DONE rev2_s2_base ${P}"
