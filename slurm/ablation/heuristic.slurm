#!/bin/bash
#SBATCH --job-name=ablation_heuristic
#SBATCH --output=logs/ablation_heuristic_%A_%a.out
#SBATCH --error=logs/ablation_heuristic_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL

set -euo pipefail

TAXON="${1:-}"
usage() {
    echo "Usage: sbatch $0 <taxon_id>" >&2
    echo "  taxon_id must be one of: 272844 (small) 759272 (medium) 10116 (large)" >&2
}
case "$TAXON" in
    272844|759272|10116) ;;
    "")
        usage
        exit 1
        ;;
    *)
        echo "Invalid taxon_id: $TAXON" >&2
        usage
        exit 1
        ;;
esac

[ -f config.env ] && source config.env

module purge
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}

mkdir -p logs

PRED_DIR=${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions
COMPLEXES=complex_coherence/protein_complexes.tsv
HIERARCHY=data/go_hierarchy_cc.tsv

# --- Dry-run (counts only, seconds to run) ---
echo "===== DRY RUN (taxon $TAXON) ====="
python -u complex_coherence/heuristic_ablation.py --dry_run \
    --predictions "$PRED_DIR/predictions_fold_01_taxon_${TAXON}.tsv" \
    --complexes "$COMPLEXES" --go_hierarchy "$HIERARCHY" \
    --taxon_id "$TAXON" --threshold 0.3 --top_k 50

# --- Full solve (with all 7 configs) ---
echo ""
echo "===== FULL SOLVE (taxon $TAXON) ====="
python -u complex_coherence/heuristic_ablation.py \
    --predictions "$PRED_DIR/predictions_fold_01_taxon_${TAXON}.tsv" \
    --complexes "$COMPLEXES" --go_hierarchy "$HIERARCHY" \
    --taxon_id "$TAXON" --threshold 0.3 --top_k 50 --timeout 600
