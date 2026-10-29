#!/bin/bash
set -euo pipefail
# Required: EXTENSION_CODE, REVISION_RESULTS, RESOURCE_DIR, EXTENSION_OUT,
# ID_MAP, COMPLEXPORTAL_DIR, PYTHON, TASKS, WORKERS. Slurm supplies worker index.
worker=${SLURM_ARRAY_TASK_ID:-0}
index=0
while read -r mode predictor taxon; do
  if (( index % WORKERS == worker )); then
    target="$EXTENSION_OUT/$mode/$predictor/$taxon.json"
    if [[ -e "$target" ]]; then
      echo "Refusing to overwrite existing result: $target" >&2
      exit 1
    fi
    "$PYTHON" "$EXTENSION_CODE/analysis/revision_extensions.py" "$mode" \
      --resources "$RESOURCE_DIR" --results "$REVISION_RESULTS" \
      --predictor "$predictor" --taxon "$taxon" --output "$target" \
      --id-map "$ID_MAP" --complexportal "$COMPLEXPORTAL_DIR"
  fi
  index=$((index+1))
done < "$TASKS"
