#!/usr/bin/env bash
# Rerun adjustment/evaluation on fixed, identified predictor inputs.
set -euo pipefail
config=$(realpath "${1:?Usage: bash reproducibility/run_all_slurm.sh /absolute/path/to/revision.env}")
source "$config"
: "${CODE_DIR:?}" "${DATA_DIR:?}" "${RESOURCE_DIR:?}" "${OUT_DIR:?}" "${PYTHON:?}" "${GAEF_DIR:?}"
: "${FUNTAXIS_MANIFEST:?}" "${PROTEOME_ROOT:?}" "${PROTEOME_INDEX:?}" "${CP_SUBUNIT_INDEX:?}"
: "${ANNOTATION_DIR:?Set the rebuilt source-reference directory}"
: "${SLURM_ACCOUNT:?}" "${SLURM_PARTITION:?}"
AUDIT_WORKERS=${AUDIT_WORKERS:-36}
AUDIT_TIME=${AUDIT_TIME:-12:00:00}
AUDIT_MEMORY=${AUDIT_MEMORY:-64G}
AUDIT_AGGREGATE_MEMORY=${AUDIT_AGGREGATE_MEMORY:-128G}
AUDIT_LINEAGE_MEMORY=${AUDIT_LINEAGE_MEMORY:-128G}
TASKS="$CODE_DIR/reproducibility/revision_tasks.txt"
[[ "$AUDIT_WORKERS" =~ ^[1-9][0-9]*$ ]]
for file in "$TASKS" "$FUNTAXIS_MANIFEST" "$PROTEOME_INDEX" "$CP_SUBUNIT_INDEX" "$ANNOTATION_DIR/reference-manifest.json"; do
  [[ -s "$file" ]] || { printf 'Missing required input: %s\n' "$file" >&2; exit 1; }
done
if [[ -e "$OUT_DIR/run-config.env" ]]; then
  printf 'Use a new output directory; a run configuration already exists.\n' >&2
  exit 1
fi
mkdir -p "$OUT_DIR/logs" "$OUT_DIR/summaries"
# Freeze configuration values. Source code and external inputs must remain fixed.
for key in CODE_DIR DATA_DIR RESOURCE_DIR OUT_DIR PYTHON GAEF_DIR FUNTAXIS_MANIFEST PROTEOME_ROOT PROTEOME_INDEX CP_SUBUNIT_INDEX ANNOTATION_DIR TASKS AUDIT_WORKERS AUDIT_MEMORY AUDIT_AGGREGATE_MEMORY AUDIT_LINEAGE_MEMORY; do
  printf '%s=%q\n' "$key" "${!key}"
done > "$OUT_DIR/run-config.env"
for model in deepgo-se seq_sim_timeset sprof-go seq_sim_heldout mlp_heldout; do
  for aspect in cc mf bp; do printf '%s %s\n' "$model" "$aspect"; done
done > "$OUT_DIR/aggregate-tasks.txt"
export AUDIT_CONFIG="$OUT_DIR/run-config.env"
common=(--parsable --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" --cpus-per-task=1 --export=ALL,AUDIT_CONFIG="$AUDIT_CONFIG")
phase="$CODE_DIR/reproducibility/revision_phase.sh"
workers=$(sbatch "${common[@]}" --array="0-$((AUDIT_WORKERS-1))" --time="$AUDIT_TIME" --mem="$AUDIT_MEMORY" --output="$OUT_DIR/logs/workers_%A_%a.log" "$phase" workers)
printf 'workers %s\n' "$workers" | tee -a "$OUT_DIR/jobs.txt"
evidence=$(sbatch "${common[@]}" --time=02:00:00 --mem=16G --output="$OUT_DIR/logs/evidence_%j.log" "$phase" evidence)
printf 'evidence %s\n' "$evidence" | tee -a "$OUT_DIR/jobs.txt"
aggregate=$(sbatch "${common[@]}" --array=0-14%4 --time=04:00:00 --mem="$AUDIT_AGGREGATE_MEMORY" --dependency="afterok:$workers" --output="$OUT_DIR/logs/aggregate_%A_%a.log" "$phase" aggregate)
printf 'aggregate %s\n' "$aggregate" | tee -a "$OUT_DIR/jobs.txt"
lineage=$(sbatch "${common[@]}" --time=04:00:00 --mem="$AUDIT_LINEAGE_MEMORY" --dependency="afterok:$workers" --output="$OUT_DIR/logs/lineage_%j.log" "$phase" lineage)
printf 'lineage %s\n' "$lineage" | tee -a "$OUT_DIR/jobs.txt"
membership=$(sbatch "${common[@]}" --array=0-5 --time=02:00:00 --mem=16G --dependency="afterok:$workers:$evidence" --output="$OUT_DIR/logs/membership_%A_%a.log" "$phase" membership)
printf 'membership %s\n' "$membership" | tee -a "$OUT_DIR/jobs.txt"
validation=$(sbatch "${common[@]}" --time=02:00:00 --mem=16G --dependency="afterok:$workers" --output="$OUT_DIR/logs/export-qc_%j.log" "$phase" validate-exports)
printf 'validate-exports %s\n' "$validation" | tee -a "$OUT_DIR/jobs.txt"
generate=$(sbatch "${common[@]}" --time=00:30:00 --mem=4G --dependency="afterok:$aggregate:$lineage:$membership:$evidence:$validation" --output="$OUT_DIR/logs/generate_%j.log" "$phase" generate)
printf 'generate %s\n' "$generate" | tee -a "$OUT_DIR/jobs.txt"
