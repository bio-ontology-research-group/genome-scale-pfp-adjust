#!/bin/bash
#SBATCH --job-name=todo-runtime-terms
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%j.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
mkdir -p $OUT
$PY -u $BASE/todo-reruns-20261031/runtime_termcentric.py --out-dir $ACCEPTED --code $CODE --go $RES/go-basic.obo --output $OUT/runtime-termcentric.json
echo DONE
