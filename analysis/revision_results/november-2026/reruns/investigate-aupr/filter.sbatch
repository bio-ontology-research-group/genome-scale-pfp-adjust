#!/bin/bash
#SBATCH --job-name=inv-aupr-filter
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-14
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
D=$BASE/todo-reruns-20261031/investigate-aupr
O=$BASE/out/investigate-aupr
mkdir -p $O
read -r model aspect < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$ACCEPTED/aggregate-tasks.txt")
cd $CODE
$PY -u $D/filter_eval.py --code $CODE --tasks $CODE/reproducibility/revision_tasks.txt --results $ACCEPTED \
  --go $RES/go-basic.obo --annotations $REFS \
  --proteins $CODE/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv \
  --heldout-train $CODE/splits/heldout/train_organisms.txt --resources $RES \
  --predictor $model --aspect $aspect --summary $ACCEPTED/summaries/${model}_${aspect}.json \
  --output $O/filter_${model}_${aspect}.json
echo DONE
