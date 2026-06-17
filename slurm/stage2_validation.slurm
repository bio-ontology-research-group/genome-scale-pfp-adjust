#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=stage2_val
#SBATCH --output=%x.log
#
# Raise curated coverage of the Stage-2 validation: run Stage-1+2 (complex
# coherence) on the organisms that have EBI Complex Portal references (E. coli,
# S. pombe, C. elegans, D. melanogaster, D. rerio, X. laevis) for DeepGO-SE and
# SPROF-GO CC, then re-validate promotions/demotions against ComplexPortal.
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

OUTROOT=${OUT_DIR}/stage2_validation
mkdir -p ${OUTROOT}
TAXA=${OUTROOT}/cp_taxa.txt
printf '83333\n284812\n6239\n7227\n7955\n8355\n' > ${TAXA}
IDX=${OUT_DIR}/../refs/complex_validation/subunit_index.json   # fallback
IDX=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/refs/complex_validation/subunit_index.json
D=${DATA_DIR}/swissprot_proteomes_folds

run() {  # predictor  tau
  local P=$1 TAU=$2
  local PRED=${D}/${P}_results/cc/predictions
  local S1=${D}/${P}_results/cc/optimized_without_taxons        # Stage-1 only (existing)
  local S2=${OUTROOT}/${P}/cc/stage2                            # Stage-1+2 (new, these taxa)
  mkdir -p ${S2}
  echo "=== ${P}: Stage-1+2 on Complex Portal taxa (tau=${TAU}) ==="
  python pipeline/run_adjustment_pipeline.py \
      --predictions_dir ${PRED} --optimized_dir ${S2} \
      --taxon_ids_file ${TAXA} \
      --optimized --complex_coherence --top_k 50 --default_threshold ${TAU} 2>&1 | tail -3
  echo "=== ${P}: validate against Complex Portal ==="
  python analysis/complex_validation/validate_stage2.py \
      --stage1_dir ${S1} --stage2_dir ${S2} \
      --subunit_index ${IDX} --threshold ${TAU} --label ${P} \
      --output_file ${OUTROOT}/${P}_cp.json 2>&1 | grep -E 'total|coverage|precision|good_demotion|harmful|n_organism' | head -20
}

run deepgo-se 0.3
run sprof-go 0.001
echo "DONE stage2_validation"
