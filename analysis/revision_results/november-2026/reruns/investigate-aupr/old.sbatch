#!/bin/bash
#SBATCH --job-name=inv-aupr-old
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=2
#SBATCH --array=0-3
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
OLD=$BASE/code
ROOT=/ibex/scratch/projects/c2014/asaad/swissprot_proteomes_folds
D=$BASE/todo-reruns-20261031/investigate-aupr
O=$BASE/out/investigate-aupr
mkdir -p $O
cd $OLD
export PYTHONPATH=$OLD:$BASE/GAEF:$OLD/pipeline
case $SLURM_ARRAY_TASK_ID in
 0) $PY -u $D/old_curve.py --oldcode $OLD --gaef $BASE/GAEF --go $OLD/data/go-basic.obo --root $ROOT --pred deepgo-se --aspect cc --output $O/old_curve_deepgo-se_cc_all.json ;;
 1) $PY -u $D/old_curve.py --oldcode $OLD --gaef $BASE/GAEF --go $OLD/data/go-basic.obo --root $ROOT --pred deepgo-se --aspect cc --taxa $D/timeset_50_taxons.txt --output $O/old_curve_deepgo-se_cc_t50.json ;;
 2) $PY -u $D/old_export_diff.py --root $ROOT --pred deepgo-se --aspect cc --tau 0.3 --output $O/old_export_diff_deepgo-se_cc.json ;;
 3) $PY -u $D/old_curve.py --oldcode $OLD --gaef $BASE/GAEF --go $OLD/data/go-basic.obo --root $ROOT --pred sprof-go --aspect cc --output $O/old_curve_sprof-go_cc_all.json ;;
esac
echo DONE
