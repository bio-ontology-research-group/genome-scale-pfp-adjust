#!/bin/bash
#SBATCH --job-name=il-termstats
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=6:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-99
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-lineage/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
D=$BASE/todo-reruns-20261031/investigate-lineage
O=$BASE/out/investigate-lineage/termstats
i=0
while read -r pred asp t tau f; do
  if (( i % 100 == SLURM_ARRAY_TASK_ID )); then
    mkdir -p $O/$pred/$asp
    out=$O/$pred/$asp/$t.json.gz
    [ -s $out ] || $PY -u $D/termstats.py --code $CODE --resources $RES --predictions $f --aspect $asp --tau $tau --output $out
  fi
  i=$((i+1))
done < $D/termstats-tasks.txt
echo DONE
