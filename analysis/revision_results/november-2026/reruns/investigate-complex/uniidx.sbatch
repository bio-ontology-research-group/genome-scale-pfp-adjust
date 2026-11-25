#!/bin/bash
#SBATCH --job-name=ic-uniidx
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-complex/logs/%x_%j.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/investigate-complex/env.sh
$PY -u $IC/build_uniprot_complex_index.py --code $CODE --go $RES/go-basic.obo --evidence $ACCEPTED/evidence-sensitivity.json --out $OI/uniprot_index
