#!/bin/bash
# Local fallback: run lineage_variants.py for every taxon of PRED whose 3 termstats exist.
# Usage: local_variants.sh PRED VARIANTS OUTSUB NPAR
set -u
W=$HOME/tmp/il; C=$HOME/Documents/papers/asaad-genomescale-paper/genome-scale-pfp-adjust
S=$HOME/Documents/papers/asaad-genomescale-paper/todo-reruns-2026-10-31/investigate-lineage
PRED=$1; VARS=$2; SUB=$3; NPAR=${4:-6}
awk -v p=$PRED '$1==p && $2=="cc"{print $3}' $W/termstats-tasks.txt | while read t; do
  [ -e $W/$SUB/.done_${PRED}_$t ] && continue
  [ -s $W/ts/$PRED/cc/$t.json.gz ] && [ -s $W/ts/$PRED/mf/$t.json.gz ] && [ -s $W/ts/$PRED/bp/$t.json.gz ] && echo $t
done | xargs -P $NPAR -I{} sh -c "PYTHONHASHSEED=0 nice -n 10 $HOME/miniconda/bin/python $S/lineage_variants.py --code $C --resources $W/res --termstats $W/ts --predictor $PRED --taxa {} --variants $VARS --truth --ic-file $W/res/go_ic.tsv --time-limit 120 --outdir $W/$SUB 2>&1 | grep -v -E '^\[INFO\]|^Computed' > $W/llogs/var_${SUB}_${PRED}_{}.log && touch $W/$SUB/.done_${PRED}_{}"
