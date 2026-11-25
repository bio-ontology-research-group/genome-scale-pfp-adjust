#!/bin/bash
#SBATCH --job-name=il-variants
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-lineage/logs/%x_%A_%a.log
# Usage: sbatch --array=0-(N-1) variants.sbatch VARIANTS OUTSUBDIR NCHUNK PRED_REGEX
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
D=$BASE/todo-reruns-20261031/investigate-lineage
O=$BASE/out/investigate-lineage
VARIANTS=$1; SUB=$2; N=${3:-100}; RE=${4:-.}
# work units: predictor taxon (588 lines)
awk -v re="$RE" '$2=="cc" && $1 ~ re {print $1, $3}' $D/termstats-tasks.txt > $O/units.txt.$SLURM_ARRAY_TASK_ID
declare -A TAXA
i=0
while read -r pred t; do
  if (( i % N == SLURM_ARRAY_TASK_ID )); then TAXA[$pred]="${TAXA[$pred]:+${TAXA[$pred]},}$t"; fi
  i=$((i+1))
done < $O/units.txt.$SLURM_ARRAY_TASK_ID
rm -f $O/units.txt.$SLURM_ARRAY_TASK_ID
for pred in "${!TAXA[@]}"; do
  $PY -u $D/lineage_variants.py --code $CODE --resources $RES --termstats $O/termstats --predictor $pred \
    --taxa "${TAXA[$pred]}" --variants "$VARIANTS" --truth \
    --ic-file $BASE/out/cost_ablation/go_ic.tsv --outdir $O/$SUB
done
echo DONE
