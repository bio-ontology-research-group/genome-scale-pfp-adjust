#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=ablation_csv
#SBATCH --output=%x_%j.log
#
# Materialize the Stage-2 scalability ablation (paper Section 3.6) to a committed
# CSV: dry-run variable/constraint counts for every heuristic configuration on the
# three benchmark proteomes (272844 small, 759272 medium, 10116 large), plus the
# real +H1+H2+H3 solve time on the large proteome. The dry-run preprocessing is
# expensive on the 20k-protein proteome, hence a batch job rather than srun.
set -euo pipefail
[ -f config.env ] && source config.env
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
export PYTHONPATH="$(pwd):${PYTHONPATH:-}"

PD=${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions
python analysis/ablation_var_counts.py \
    --predictions_dir ${PD} \
    --output analysis/revision_results/ablation_var_counts.csv \
    --solve_taxon 10116
echo "DONE ablation_var_counts"
