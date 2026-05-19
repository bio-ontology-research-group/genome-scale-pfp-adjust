#!/usr/bin/env bash
# run_term_flip_ranking.sh
#
# Run term-centric flip ranking for all 9 method x ontology combinations,
# both taxon-guided (optimized/) and taxon-unknown (optimized_without_taxons/)
# conditions, using the timeset test proteins.
#
# A "flip" is counted only when the ORIGINAL prediction score exceeds
# $SCORE_THRESHOLD (default 0.1). Flips below that are dismissed as
# low-confidence noise.
#
# Prerequisites:
#   conda activate ${CONDA_ENV}
#   Run from: repository root
#
# Outputs land in: term-flip-ranking/
#   <method>_<ont>_<cond>_by_count.tsv   -- top-K terms by flip count
#   <method>_<ont>_<cond>_by_cost.tsv    -- top-K terms by cumulative cost
#   <method>_<ont>_<cond>_full.tsv       -- every term with >=1 flip
#   <method>_<ont>_<cond>_summary.tsv    -- totals
#   run_log.txt

set -euo pipefail

SCRIPT="analysis/term_flip_ranking.py"
TEST_PROTEINS="splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv"
ANNOTATIONS_DIR="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic"
GO_FILE="data/go-basic.obo"
TAXON_CONSTRAINTS="data/go_taxon_constraints_extracted_obo.tsv"
OUT_BASE="term-flip-ranking"
LOG="${OUT_BASE}/run_log.txt"
SCORE_THRESHOLD="${SCORE_THRESHOLD:-0.1}"
TOP_K="${TOP_K:-25}"

mkdir -p "${OUT_BASE}"
: > "${LOG}"

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
            PREFIX="${OUT_BASE}/${TAG}"

            echo "=== ${TAG} (threshold=${SCORE_THRESHOLD}) ===" | tee -a "${LOG}"
            python3 "${SCRIPT}" \
                --predictions_dir      "${ROOT}/${ont}/predictions" \
                --optimized_dir        "${OPT_DIR}" \
                --test_proteins_file   "${TEST_PROTEINS}" \
                --go_file              "${GO_FILE}" \
                --subontology          "${ont}" \
                --output_prefix        "${PREFIX}" \
                --score_threshold      "${SCORE_THRESHOLD}" \
                --top_k                "${TOP_K}" \
                --taxon_constraints_file "${TAXON_CONSTRAINTS}" \
                --annotations_dir      "${ANNOTATIONS_DIR}" \
                2>&1 | tee -a "${LOG}"
            echo "" | tee -a "${LOG}"
        done
    done
done

echo "Done. Results in ${OUT_BASE}/"
