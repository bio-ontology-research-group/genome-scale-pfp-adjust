#!/bin/bash
#SBATCH --job-name=todo-contigs
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=48:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-14
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
D=$BASE/todo-reruns-20261031/contigs
DATA=$BASE/baseline-rebuild-20261029/predictions
TAXA=(83332 3572 69014)
COMBOS=("seq_sim_timeset bp" "deepgo-se bp" "sprof-go bp" "seq_sim_timeset cc" "seq_sim_timeset mf")
T=${TAXA[$((SLURM_ARRAY_TASK_ID / 5))]}
read PRED ASP <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID % 5))]}"
TAU=0.001; [ "$PRED" = deepgo-se ] && TAU=0.3
SIZE=$(grep 'Total concatenated' $D/cds_$T.log | sed -E 's/.*size: ([0-9,]+) bp.*/\1/' | tr -d ,)
F=$(ls $DATA/${PRED}_results/$ASP/predictions/predictions_fold_*_taxon_${T}.tsv)
mkdir -p $D/results
$PY -u $BASE/todo-reruns-20261031/contig_lineage.py --code $CODE --resources $RES --predictions $F \
  --cds $D/cds_$T.tsv --genome-size $SIZE --taxon $T --aspect $ASP --tau $TAU \
  --output $D/results/${T}_${PRED}_${ASP}.json
echo DONE
