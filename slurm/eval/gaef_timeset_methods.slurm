#!/bin/bash
#SBATCH --job-name=eval_gaef_timeset
#SBATCH --output=logs/eval_gaef_timeset_%A_%a.out
#SBATCH --error=logs/eval_gaef_timeset_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --array=0-8
#
# GAEF evaluation (evaluate_directory.py) for timeset predictions, array layout:
#   9 tasks = 3 methods × 3 ontologies (cc, mf, bp).
#   Each task runs:
#     - .../<ont>/predictions           -> eval_gaef_predictions
#     - .../<ont>/optimized_without_taxons -> eval_gaef_optimized_without_taxons
#   For ont=cc only (same task): if present,
#     - .../cc/optimized_without_taxons_complex_coherence -> eval_gaef_optimized_without_taxons_complex_coherence
#
# Index: TASK_ID = method_index * 3 + ontology_index
#   method 0: seq-sim     ontology 0,1,2 -> cc, mf, bp
#   method 1: DeepGO-SE
#   method 2: SPROF-GO
#
# Submit from the repository root so relative paths resolve:
#   sbatch slurm/eval/gaef_timeset_methods.slurm
#
# Optional env overrides (same shell):
#   BASE_DIR=/path/to/swissprot_proteomes_folds sbatch slurm/eval/gaef_timeset_methods.slurm
#
# Run a subset of array indices (e.g. only seq-sim):
#   sbatch --array=0-2 slurm/eval/gaef_timeset_methods.slurm

set -euo pipefail

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

BASE_DIR="${BASE_DIR:-${DATA_DIR}/swissprot_proteomes_folds}"
ANN_DIR="${BASE_DIR}/annotations-go-basic"
GO_FILE="${GO_OBO:-data/go-basic.obo}"
GAEF_DIR="${GAEF_DIR:?GAEF_DIR not set; source config.env}"

SEQ_SIM_ROOT="${SEQ_SIM_ROOT:-${BASE_DIR}/seq_sim_timeset_results}"
DEEPGO_TIMESET_ROOT="${DEEPGO_TIMESET_ROOT:-${BASE_DIR}/deepgo-se_results}"
SPROF_TIMESET_ROOT="${SPROF_TIMESET_ROOT:-${BASE_DIR}/sprof-go_results}"

THR_SEQ_SIM="${THR_SEQ_SIM:-0.3}"
THR_DEEPGO="${THR_DEEPGO:-0.3}"
THR_SPROF="${THR_SPROF:-0.3}"

threshold_args() {
  local ROOT_SUB="$1"
  local DEF_THR="$2"
  if [[ -f "${ROOT_SUB}/fold_thresholds.json" ]]; then
    echo --threshold_file "${ROOT_SUB}/fold_thresholds.json"
  else
    echo --threshold "${DEF_THR}"
  fi
}

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID not set — run under sbatch --array=0-8}"
METHOD_IDX=$((TASK_ID / 3))
ONTO_IDX=$((TASK_ID % 3))
ONTOLOGIES=(cc mf bp)
TASK="${ONTOLOGIES[$ONTO_IDX]}"

case "${METHOD_IDX}" in
  0)
    NAME="seq-sim (timeset)"
    MROOT="${SEQ_SIM_ROOT}"
    DEF_THR="${THR_SEQ_SIM}"
    ;;
  1)
    NAME="DeepGO-SE (timeset)"
    MROOT="${DEEPGO_TIMESET_ROOT}"
    DEF_THR="${THR_DEEPGO}"
    ;;
  2)
    NAME="SPROF-GO (timeset)"
    MROOT="${SPROF_TIMESET_ROOT}"
    DEF_THR="${THR_SPROF}"
    ;;
  *)
    echo "Invalid METHOD_IDX=${METHOD_IDX} from TASK_ID=${TASK_ID}"
    exit 1
    ;;
esac

echo "========================================"
echo "Array task ${TASK_ID} / job ${SLURM_ARRAY_JOB_ID:-local}"
echo "Method [${METHOD_IDX}]: ${NAME}"
echo "Ontology [${ONTO_IDX}]: ${TASK}"
echo "Root: ${MROOT}"
echo "========================================"

if [[ ! -d "${MROOT}" ]]; then
  echo "[Error] Method root not found: ${MROOT}"
  exit 1
fi

TA="${MROOT}/${TASK}"
TARGS=$(threshold_args "${TA}" "${DEF_THR}")

echo "--- ${NAME} ${TASK}: predictions (before) ---"
python pipeline/evaluate_directory.py \
  --predictions_dir "${TA}/predictions" \
  --annotations_dir "${ANN_DIR}" \
  --output_dir "${TA}/eval_gaef_predictions" \
  ${TARGS} \
  --go_file "${GO_FILE}" \
  --GAEF_dir "${GAEF_DIR}" \
  --subontology "${TASK}"

echo "--- ${NAME} ${TASK}: optimized_without_taxons ---"
python pipeline/evaluate_directory.py \
  --predictions_dir "${TA}/optimized_without_taxons" \
  --annotations_dir "${ANN_DIR}" \
  --output_dir "${TA}/eval_gaef_optimized_without_taxons" \
  ${TARGS} \
  --go_file "${GO_FILE}" \
  --GAEF_dir "${GAEF_DIR}" \
  --subontology "${TASK}"

if [[ "${TASK}" == "cc" ]]; then
  CC_COH_DIR="${MROOT}/cc/optimized_without_taxons_complex_coherence"
  if [[ -d "${CC_COH_DIR}" ]]; then
    echo "--- ${NAME} cc: optimized_without_taxons_complex_coherence ---"
    TARGS_CC=$(threshold_args "${MROOT}/cc" "${DEF_THR}")
    python pipeline/evaluate_directory.py \
      --predictions_dir "${CC_COH_DIR}" \
      --annotations_dir "${ANN_DIR}" \
      --output_dir "${MROOT}/cc/eval_gaef_optimized_without_taxons_complex_coherence" \
      ${TARGS_CC} \
      --go_file "${GO_FILE}" \
      --GAEF_dir "${GAEF_DIR}" \
      --subontology cc
  else
    echo "[Skip] ${NAME}: missing ${CC_COH_DIR}"
  fi
fi

echo "Done. Summaries under ${TA}/eval_gaef_*/evaluation_summary.json (and cc complex-coherence when run)."
