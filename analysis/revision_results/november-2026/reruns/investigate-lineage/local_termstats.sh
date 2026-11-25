#!/bin/bash
# Local fallback when the IBEX queue is saturated: fetch one prediction file
# (rsync -z), summarize it with termstats.py, delete it. Usage: local_termstats.sh PRED WORKER NWORKERS
set -u
W=$HOME/tmp/il; C=$HOME/Documents/papers/asaad-genomescale-paper/genome-scale-pfp-adjust
S=$HOME/Documents/papers/asaad-genomescale-paper/todo-reruns-2026-10-31/investigate-lineage
PRED=$1; K=$2; N=$3; i=0
mkdir -p $W/fetch/$K
awk -v p=$PRED '$1==p' $W/termstats-tasks.txt | while read -r pred asp t tau f; do
  if (( i % N == K )); then
    out=$W/ts/$pred/$asp/$t.json.gz; mkdir -p $(dirname $out)
    if [ ! -s $out ]; then
      loc=$W/pred/$pred/$asp/$(basename $f)
      if [ ! -s $loc ] || [ $(stat -c %s $loc) -lt 1000 ]; then
        loc=$W/fetch/$K/$(basename $f); rsync -az ibex:$f $loc 2>/dev/null
      fi
      PYTHONHASHSEED=0 nice -n 10 $HOME/miniconda/bin/python $S/termstats.py --code $C --resources $W/res --predictions $loc --aspect $asp --tau $tau --output $out > /dev/null && echo "$pred $asp $t ok"
      [[ $loc == $W/fetch/* ]] && rm -f $loc
    fi
  fi
  i=$((i+1))
done
