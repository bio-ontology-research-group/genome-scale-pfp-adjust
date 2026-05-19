#!/bin/bash
#SBATCH --job-name=uniprot_predict_sprof_go
#SBATCH --output=logs/uniprot_predict_sprof_go_%j.out
#SBATCH --error=logs/uniprot_predict_sprof_go_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --constraint=[v100|rtx2080ti|gtx1080ti]
#SBATCH --mail-type=ALL

# Runs the pretrained SPROF-GO model on the timeset test organisms, then converts
# its raw outputs to the per-protein TSV format consumed by the solver.
# Writes ${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results/{cc,mf,bp}/predictions/.

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

WORK_DIR=${WORK_DIR:-${DATA_DIR}/sprof_go/sprof_go_work}
OUTPUT_DIR=${OUTPUT_DIR:-${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results}
mkdir -p ${WORK_DIR} ${OUTPUT_DIR}

# SPROF-GO inference command template. Override with SPROF_CMD if your install differs.
# The {fasta} and {output_dir} placeholders are filled in by run_sprof_go.py.
SPROF_CMD="${SPROF_CMD:-python ${SPROF_GO_REPO}/predict.py --input_fasta {fasta} --output_dir {output_dir}}"

# Stage 1: run SPROF-GO for each timeset test organism
python baselines/sprof_go/run_sprof_go.py \
    --test-organisms-file splits/timeset/test_organisms.txt \
    --proteomes-dir ${DATA_DIR}/uniprot_reference_proteomes \
    --proteomes-ids data/uniprot_proteomes_ids.tsv \
    --output-dir ${OUTPUT_DIR} \
    --work-dir ${WORK_DIR} \
    --subontology all \
    --sprof-cmd "${SPROF_CMD}"

# Stage 2: convert raw SPROF-GO outputs to the canonical prediction TSV layout
python baselines/sprof_go/process_sprof-go_results.py \
    --work-dir ${WORK_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --subontology all
