#!/bin/bash
# Usage: aggregate_full.sh TASKFILE NPAR PRED...   (waits for commit headroom NEED_GB, retries up to 4x)
source /mnt/data1/leechuck/genoadjust-ws/ws.env
export CODE PYTHONPATH=$PYTHONPATH:$W/rerun-notaxon PYTHONHASHSEED
R=$W/rerun-notaxon; T=$1; NPAR=$2; shift 2; mkdir -p $R/summaries_full $R/logs/agg_full
export R T PY RES REFS NEED_GB=${NEED_GB:-24}
for p in "$@"; do for a in cc mf bp; do echo $p $a; done; done | xargs -P $NPAR -L 1 sh -c '
  for try in 1 2 3 4; do
    while [ $(awk "/CommitLimit/{l=\$2}/Committed_AS/{c=\$2}END{print int((l-c)/1048576)}" /proc/meminfo) -lt $NEED_GB ]; do sleep $((30 + $$ % 30)); done
    if nice -n 10 $PY $R/scripts/notaxon_aggregate_full.py --tasks $R/$T --results $R/out_full --go $RES/go-basic.obo --annotations $REFS \
     --proteins $CODE/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv --heldout-train $CODE/splits/heldout/train_organisms.txt \
     --predictor $0 --aspect $1 --output $R/summaries_full/$0_$1.json > $R/logs/agg_full/$0_$1.log 2>&1; then echo "AGG OK $0 $1"; break
    else echo "AGG FAIL $0 $1 try $try"; sleep 90; fi
  done'
echo AGG DONE
