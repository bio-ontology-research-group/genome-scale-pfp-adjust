#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=genoadjust-audit
set -euo pipefail
: "${AUDIT_CONFIG:?Set AUDIT_CONFIG to the experiment configuration file}"
source "$AUDIT_CONFIG"
: "${CODE_DIR:?}" "${DATA_DIR:?}" "${RESOURCE_DIR:?}" "${OUT_DIR:?}" "${PYTHON:?}" "${TASKS:?}"
cd "$CODE_DIR"
export PYTHONPATH="$CODE_DIR:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1
export TMPDIR="$OUT_DIR/scratch/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"
read -r predictor aspect taxon < <(sed -n "$(( ${SLURM_ARRAY_TASK_ID:-0} + 1 ))p" "$TASKS")
[[ -n "$taxon" ]]
options=(--infer-taxons)
if [[ -n "${ANNOTATION_DIR:-}" ]]; then options+=(--annotations "$ANNOTATION_DIR"); fi
if [[ "$predictor" != seq_sim_heldout ]]; then options+=(--sweep); fi
if [[ "$predictor" != *_heldout && "$aspect" == cc ]]; then options+=(--stage2 --soft); fi
if [[ "$predictor" != *_heldout ]]; then
    options+=(--proteins "$CODE_DIR/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv")
fi
if [[ -n "${FUNTAXIS_MANIFEST:-}" && "$predictor" != *_heldout ]]; then
    options+=(--funtaxis "$FUNTAXIS_MANIFEST")
fi
command=("$PYTHON" analysis/revision_audit_run.py \
    --data "$DATA_DIR" --resources "$RESOURCE_DIR" \
    --predictor "$predictor" --aspect "$aspect" --taxon "$taxon" \
    --output "$OUT_DIR/$predictor/$aspect/$taxon.json.gz" "${options[@]}")
attempt_log="$OUT_DIR/logs/attempt_${predictor}_${aspect}_${taxon}_120.log"
mkdir -p "$OUT_DIR/logs"
if "${command[@]}" --solver-seconds 120 2>&1 | tee "$attempt_log"; then
    exit 0
fi
# Retry only a solver timeout without a feasible incumbent. Every other failure
# remains fatal. The accepted JSON records the actual budget and all attempts
# remain in the logs; no partial result is accepted.
if [[ "$aspect" == cc && "$predictor" != *_heldout ]] && \
   grep -q '^RuntimeError: CP-SAT did not produce a valid repair: UNKNOWN$' "$attempt_log"; then
    "${command[@]}" --solver-seconds 900 2>&1 | tee "$OUT_DIR/logs/attempt_${predictor}_${aspect}_${taxon}_900.log"
else
    exit 1
fi
