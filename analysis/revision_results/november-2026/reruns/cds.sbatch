#!/bin/bash
#SBATCH --job-name=todo-cds
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=3
#SBATCH --output=/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/logs/%x_%j.log
set -euo pipefail
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
P=/ibex/scratch/projects/c2014/asaad/uniprot_reference_proteomes
D=$BASE/todo-reruns-20261031/contigs
cd $D
run() { $PY -u $CODE/analysis/metagenomics_find_locations.py --dat-file $2 --genome-file $3 --output-file $D/cds_$1.tsv > $D/cds_$1.log 2>&1; grep -E 'Total concatenated|Rows written|not in' $D/cds_$1.log; }
run 83332 $P/Bacteria/UP000001584/UP000001584_83332.dat.gz $D/AL123456.3.embl &
run 69014 $P/Archaea/UP000000536/UP000000536_69014.dat.gz $P/genomes_embl/AP006878.1.embl &
run 3572 $P/Eukaryota/UP001443914/UP001443914_3572.dat.gz $P/genomes_embl/GCA_040167595.1.embl &
wait
echo DONE
