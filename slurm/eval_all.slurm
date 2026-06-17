#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=eval_all
#SBATCH --array=0-11
#SBATCH --output=%x_%a.log
#
# Per-organism CAFA eval (predictions vs taxon-guided optimized) for the cells
# not yet covered: timeset BP/MF (3 predictors) and the full heldout split
# (seq_sim_heldout, mlp_heldout) on CC/MF/BP. Feeds the BP/MF/heldout Wilcoxon.
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):${GAEF_DIR}:${PYTHONPATH:-}"

# combo list: "predictorRoot:ontology"
COMBOS=(
  deepgo-se:mf deepgo-se:bp
  seq_sim_timeset:mf seq_sim_timeset:bp
  sprof-go:mf sprof-go:bp
  seq_sim_heldout:cc seq_sim_heldout:mf seq_sim_heldout:bp
  mlp_heldout:cc mlp_heldout:mf mlp_heldout:bp
)
SPEC=${COMBOS[$SLURM_ARRAY_TASK_ID]}
P=${SPEC%%:*}; O=${SPEC##*:}
BASE=${DATA_DIR}/swissprot_proteomes_folds/${P}_results/${O}
ANN=${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic
OUT=${OUT_DIR}/eval/${P}/${O}
mkdir -p ${OUT}

for S in predictions optimized; do
    echo "=== eval ${P}/${O} ${S} ==="
    python pipeline/evaluate_directory.py \
        --predictions_dir ${BASE}/${S} \
        --annotations_dir ${ANN} \
        --output_dir ${OUT}/${S} \
        --threshold 0.3 --subontology ${O} \
        --GAEF_dir ${GAEF_DIR} --go_file ${GO_OBO} --num_workers 8 2>&1 | tail -3
done
echo "DONE ${P}/${O}"
