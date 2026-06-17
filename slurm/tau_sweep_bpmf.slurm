#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=tausweep_bpmf
#SBATCH --array=0-5
#SBATCH --output=%x_%a.log
#
# tau sweep (Stage 1, taxon-guided) flip rates for the BP and MF sub-ontologies,
# extending the CC sweep. Array index -> (predictor, ontology).
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

PREDS=(deepgo-se seq_sim_timeset sprof-go)
ONTS=(mf bp)
P=${PREDS[$((SLURM_ARRAY_TASK_ID / 2))]}
O=${ONTS[$((SLURM_ARRAY_TASK_ID % 2))]}
BASE=${DATA_DIR}/swissprot_proteomes_folds/${P}_results/${O}/predictions
OUTROOT=${OUT_DIR}/tau_sweep/${P}/${O}
RES=${OUT_DIR}/tau_sweep/${P}_${O}_fliprates.tsv
mkdir -p ${OUTROOT}
: > ${RES}

for TAU in 0.1 0.2 0.3 0.4 0.5; do
    OD=${OUTROOT}/optimized_tau${TAU}
    echo "[${P}/${O}] tau=${TAU}"
    python pipeline/run_adjustment_pipeline.py \
        --predictions_dir ${BASE} \
        --optimized_dir   ${OD} \
        --optimized --default_threshold ${TAU} --provide_taxon_id \
        --skip_existing 2>&1 | tail -2
    python analysis/flip_rate.py \
        --predictions_dir ${BASE} \
        --optimized_dir   ${OD} \
        --threshold ${TAU} --label ${P}_${O} | tee -a ${RES}
done
echo "DONE ${P}_${O}"
