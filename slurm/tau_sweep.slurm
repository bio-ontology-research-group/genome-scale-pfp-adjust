#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=tau_sweep
#SBATCH --array=0-2
#SBATCH --output=%x_%a.log
#
# tau sensitivity sweep (Stage 1, taxon-guided) on CC for the three timeset
# predictors. One array task per predictor; loops tau in {0.1,0.2,0.3,0.4,0.5}
# and records the aggregate flip rate at each tau. Answers the reviewer
# question: are flip rates qualitatively unchanged across tau?
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

PREDS=(deepgo-se seq_sim_timeset sprof-go)
P=${PREDS[$SLURM_ARRAY_TASK_ID]}
BASE=${DATA_DIR}/swissprot_proteomes_folds/${P}_results/cc/predictions
OUTROOT=${OUT_DIR}/tau_sweep/${P}
RES=${OUT_DIR}/tau_sweep/${P}_cc_fliprates.tsv
mkdir -p ${OUTROOT}
: > ${RES}

for TAU in 0.1 0.2 0.3 0.4 0.5; do
    OD=${OUTROOT}/cc_optimized_tau${TAU}
    echo "[${P}] tau=${TAU} -> adjusting"
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE} \
        --optimized_dir   ${OD} \
        --optimized --default_threshold ${TAU} --provide_taxon_id \
        --skip_existing 2>&1 | tail -3
    python analysis/flip_rate.py \
        --predictions_dir ${BASE} \
        --optimized_dir   ${OD} \
        --threshold ${TAU} --label ${P}_cc | tee -a ${RES}
done
echo "DONE ${P}"
