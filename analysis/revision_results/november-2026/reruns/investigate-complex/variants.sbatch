#!/bin/bash
#SBATCH --job-name=ic-var
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-149%60
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/investigate-complex/logs/%x_%A_%a.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/investigate-complex/env.sh
P=${PREDS[$((SLURM_ARRAY_TASK_ID / 50))]}; T=${TAXA[$((SLURM_ARRAY_TASK_ID % 50))]}
cd $IC
$PY -u variants.py --code $CODE --go $RES/go-basic.obo --accepted $ACCEPTED --predictor $P --taxon $T --out $OI/variants
if [ "$P" = deepgo-se ]; then
  $PY -u variants.py --code $CODE --go $RES/go-basic.obo --accepted $ACCEPTED --predictor $P --taxon $T --out $OI/variants_floor015 \
    --dense-floor 0.15 --variants ev_hard,ev_prefer,ev_abstain
fi
echo DONE
