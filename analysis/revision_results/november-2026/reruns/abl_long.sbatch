#!/bin/bash
#SBATCH --job-name=todo-ablation-long
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%j.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
F=$(ls $ACCEPTED/seq_sim_timeset/cc/stage1/optimized_fold_*_taxon_272844.tsv)
$PY -u $BASE/todo-reruns-20261031/stage2_ablation.py --code $CODE --go $RES/go-basic.obo --stage1 $F \
  --taxon 272844 --tau 0.3 --timeout 7200 --workers 16 --solve-naive --output $OUT/ablation/ablation_272844_tau0.3_long.json
echo DONE
