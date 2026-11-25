#!/bin/bash
#SBATCH --job-name=ic-diag
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-149%40
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-complex/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/investigate-complex/env.sh
P=${PREDS[$((SLURM_ARRAY_TASK_ID / 50))]}; T=${TAXA[$((SLURM_ARRAY_TASK_ID % 50))]}
cd $IC
$PY -u diag.py --code $CODE --go $RES/go-basic.obo --accepted $ACCEPTED --cp_index $CPIDX --uniprot_dir $OI/uniprot_index \
  --idmap $IDMAP --predictor $P --taxon $T --out $OI/diag/$P/$T.json.gz
