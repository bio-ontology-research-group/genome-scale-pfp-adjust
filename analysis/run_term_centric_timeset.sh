#!/usr/bin/env bash
# run_term_centric_timeset.sh
#
# Re-runs the term-centric evaluation for all 9 method×ontology combinations,
# both taxon-guided (optimized/) and taxon-unknown (optimized_without_taxons/)
# conditions, using the timeset test proteins.
#
# Prerequisites:
#   conda activate ${CONDA_ENV}  (or whichever env has the pipeline dependencies)
#   Run from: repository root
#
# Outputs land in: term-centric-analysis-v2/
#   <method>_<ont>_taxons.tsv          -- top-50 terms (taxon-guided)
#   <method>_<ont>_without_taxons.tsv  -- top-50 terms (taxon-unknown)
#   <method>_<ont>_taxons_full.tsv     -- full per-term detail
#   <method>_<ont>_without_taxons_full.tsv
#   <method>_<ont>_taxons_binned.tsv   -- binned summary
#   <method>_<ont>_without_taxons_binned.tsv
#   run_log.txt                        -- flip-precision summary for every run

set -euo pipefail

SCRIPT="analysis/term_metric_diff_report.py"
TEST_PROTEINS="splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv"
ANNOTATIONS_DIR="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic"
GO_FILE="data/go-basic.obo"
OUT_BASE="term-centric-analysis-v2"
LOG="${OUT_BASE}/run_log.txt"

mkdir -p "${OUT_BASE}"
> "${LOG}"

declare -A METHOD_DIRS=(
    [seq_sim]="${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results"
    [deepgo]="${DATA_DIR}/swissprot_proteomes_folds/deepgo-se_results"
    [sprof]="${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results"
)

ONTS=(cc mf bp)

for method in seq_sim deepgo sprof; do
    ROOT="${METHOD_DIRS[$method]}"
    for ont in "${ONTS[@]}"; do
        for cond in taxons without_taxons; do
            if [[ "$cond" == "taxons" ]]; then
                OPT_DIR="${ROOT}/${ont}/optimized"
            else
                OPT_DIR="${ROOT}/${ont}/optimized_without_taxons"
            fi

            if [[ ! -d "$OPT_DIR" ]]; then
                echo "SKIP (no dir): ${method} ${ont} ${cond}" | tee -a "${LOG}"
                continue
            fi

            TAG="${method}_${ont}_${cond}"
            OUT_TOP="${OUT_BASE}/${TAG}.tsv"
            OUT_DIR="${OUT_BASE}/${TAG}_detail"

            echo "=== ${TAG} ===" | tee -a "${LOG}"
            python3 "${SCRIPT}" \
                --predictions_dir  "${ROOT}/${ont}/predictions" \
                --optimized_dir    "${OPT_DIR}" \
                --test_proteins_file "${TEST_PROTEINS}" \
                --annotations_dir  "${ANNOTATIONS_DIR}" \
                --go_file          "${GO_FILE}" \
                --subontology      "${ont}" \
                --output_file      "${OUT_TOP}" \
                --output_dir       "${OUT_DIR}" \
                --top_k            50 \
                2>&1 | tee -a "${LOG}"
            echo "" | tee -a "${LOG}"
        done
    done
done

echo "Done. Results in ${OUT_BASE}/"
