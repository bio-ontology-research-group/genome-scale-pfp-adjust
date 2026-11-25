#!/bin/bash
#SBATCH --job-name=todo-gt
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%j.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
mkdir -p $OUT/gt-gaef/reports
# prediction collection taxa (564)
ls /ibex/scratch/projects/c2014/asaad/swissprot_proteomes_folds/deepgo-se_results/cc/predictions | sed -E 's/.*_taxon_([0-9]+)\.tsv/\1/' | sort -u > $OUT/gt-gaef/collection-taxa.txt
wc -l $OUT/gt-gaef/collection-taxa.txt
cd $CODE
$PY -u pipeline/evaluate_annotations_gaef.py --annotations_dir $REFS --output_dir $OUT/gt-gaef/reports \
  --go_file $RES/go-basic.obo --GAEF_dir $BASE/GAEF --num_workers 16
$PY -u $BASE/todo-reruns-20261031/gt_taxon_audit.py --code $CODE --resources $RES --annotations $REFS \
  --taxa $OUT/gt-gaef/collection-taxa.txt --output $OUT/gt-gaef/closure-audit.json
echo DONE
