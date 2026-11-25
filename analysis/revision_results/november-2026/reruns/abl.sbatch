#!/bin/bash
#SBATCH --job-name=todo-ablation
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --mem=240G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-2
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
TAXA=(272844 759272 10116)
T=${TAXA[$SLURM_ARRAY_TASK_ID]}
NAIVE=""; [ "$T" != 10116 ] && NAIVE=--solve-naive
mkdir -p $OUT/ablation
F=$(ls $ACCEPTED/seq_sim_timeset/cc/stage1/optimized_fold_*_taxon_${T}.tsv)
$PY -u $BASE/todo-reruns-20261031/stage2_ablation.py --code $CODE --go $RES/go-basic.obo --stage1 $F \
  --taxon $T --tau 0.001 --timeout 600 $NAIVE --output $OUT/ablation/ablation_${T}.json
echo DONE
