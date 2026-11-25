#!/bin/bash
#SBATCH --job-name=il-contigs
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-14
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-lineage/logs/%x_%A_%a.log
# Usage: sbatch contigs.sbatch VARIANTS OUTSUB  (same combos/sampling as ../contigs.sbatch; 100 kb and 1 Mb kept)
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
D=$BASE/todo-reruns-20261031/contigs
DATA=$BASE/baseline-rebuild-20261029/predictions
O=$BASE/out/investigate-lineage/$2; mkdir -p $O
TAXA=(83332 3572 69014)
COMBOS=("seq_sim_timeset bp" "deepgo-se bp" "sprof-go bp" "seq_sim_timeset cc" "seq_sim_timeset mf")
T=${TAXA[$((SLURM_ARRAY_TASK_ID / 5))]}
read PRED ASP <<< "${COMBOS[$((SLURM_ARRAY_TASK_ID % 5))]}"
TAU=0.001; [ "$PRED" = deepgo-se ] && TAU=0.3
SIZE=$(grep 'Total concatenated' $D/cds_$T.log | sed -E 's/.*size: ([0-9,]+) bp.*/\1/' | tr -d ,)
F=$(ls $DATA/${PRED}_results/$ASP/predictions/predictions_fold_*_taxon_${T}.tsv)
$PY -u $BASE/todo-reruns-20261031/investigate-lineage/contig_variants.py --code $CODE --resources $RES --predictions $F \
  --cds $D/cds_$T.tsv --genome-size $SIZE --taxon $T --aspect $ASP --tau $TAU --variants "$1" \
  --output $O/${T}_${PRED}_${ASP}.json
echo DONE
