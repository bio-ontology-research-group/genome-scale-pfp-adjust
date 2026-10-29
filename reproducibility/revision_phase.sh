#!/usr/bin/env bash
set -euo pipefail
: "${AUDIT_CONFIG:?}"
source "$AUDIT_CONFIG"
export GAEF_DIR
export PYTHONPATH="$CODE_DIR:${PYTHONPATH:-}"
export TMPDIR="$OUT_DIR/scratch/${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"
cd "$CODE_DIR"
phase=${1:?Pass workers, aggregate, evidence, lineage, membership or generate}
case "$phase" in
workers)
  failures=0
  total=$(wc -l < "$TASKS")
  worker=${SLURM_ARRAY_TASK_ID:-0}
  for ((index=worker;index<total;index+=AUDIT_WORKERS)); do
    if SLURM_ARRAY_TASK_ID="$index" bash "$CODE_DIR/slurm/revision_audit.slurm"; then
      printf 'DONE task %s\n' "$index"
    else
      printf 'FAILED task %s\n' "$index" >&2
      failures=$((failures+1))
    fi
  done
  [[ "$failures" -eq 0 ]]
  ;;
aggregate)
  read -r model aspect < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$OUT_DIR/aggregate-tasks.txt")
  "$PYTHON" analysis/aggregate_revision_audit.py --tasks "$TASKS" --results "$OUT_DIR" \
    --go "$RESOURCE_DIR/go-basic.obo" --annotations "$ANNOTATION_DIR" \
    --proteins "$CODE_DIR/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv" \
    --heldout-train "$CODE_DIR/splits/heldout/train_organisms.txt" --predictor "$model" --aspect "$aspect" \
    --output "$OUT_DIR/summaries/${model}_${aspect}.json"
  ;;
evidence)
  "$PYTHON" analysis/audit_annotation_evidence.py --root "$PROTEOME_ROOT" --proteomes "$PROTEOME_INDEX" \
    --taxa "$CODE_DIR/splits/timeset/test_organisms.txt" --go "$RESOURCE_DIR/go-basic.obo" \
    --mapping "$OUT_DIR/protein-id-map.json" --output "$OUT_DIR/evidence-sensitivity.json"
  ;;
lineage)
  "$PYTHON" analysis/audit_lineage_results.py --results "$OUT_DIR" --resources "$RESOURCE_DIR" \
    --taxa "$CODE_DIR/splits/timeset/test_organisms.txt" --heldout-taxa "$CODE_DIR/splits/heldout/test_organisms.txt" \
    --output "$OUT_DIR/lineage.json"
  ;;
membership)
  models=(seq_sim_timeset deepgo-se sprof-go)
  modes=(standard demotion)
  model=${models[$((SLURM_ARRAY_TASK_ID/2))]}
  mode=${modes[$((SLURM_ARRAY_TASK_ID%2))]}
  threshold=.001
  if [[ "$model" == deepgo-se ]]; then threshold=.3; fi
  "$PYTHON" analysis/complex_validation/validate_stage2.py \
    --stage1_dir "$OUT_DIR/$model/cc/stage1" --stage2_dir "$OUT_DIR/$model/cc/stage2_$mode" \
    --subunit_index "$CP_SUBUNIT_INDEX" --protein-id-map "$OUT_DIR/protein-id-map.json" \
    --go "$RESOURCE_DIR/go-basic.obo" \
    --threshold "$threshold" --label "$model $mode" --output_file "$OUT_DIR/membership_${model}_${mode}.json"
  ;;
validate-exports)
  "$PYTHON" analysis/validate_repair_exports.py --results "$OUT_DIR" --resources "$RESOURCE_DIR" \
    --taxa "$CODE_DIR/splits/timeset/test_organisms.txt" --output "$OUT_DIR/export-validation.json"
  ;;
generate)
  "$PYTHON" reproducibility/reproduce_all_tables_and_figures.py --results "$OUT_DIR" --output "$OUT_DIR/generated"
  ;;
*) printf 'Unknown phase: %s\n' "$phase" >&2; exit 2 ;;
esac
