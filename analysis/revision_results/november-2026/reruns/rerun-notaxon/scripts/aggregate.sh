#!/bin/bash
# Usage: aggregate.sh TASKFILE NPAR PRED...
source /mnt/data1/leechuck/genoadjust-ws/ws.env
export CODE PYTHONPATH=$PYTHONPATH:$W/rerun-notaxon PYTHONHASHSEED
R=$W/rerun-notaxon; T=$1; NPAR=$2; shift 2; mkdir -p $R/summaries $R/logs/agg
for p in "$@"; do for a in cc mf bp; do echo $p $a; done; done | xargs -P $NPAR -L 1 sh -c '
  nice -n 10 '$PY' '$R'/scripts/notaxon_aggregate.py --tasks '$R/$T' --results '$R'/out --go '$RES'/go-basic.obo --annotations '$REFS' \
   --proteins '$CODE'/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv --heldout-train '$CODE'/splits/heldout/train_organisms.txt \
   --predictor $0 --aspect $1 --output '$R'/summaries/$0_$1.json > '$R'/logs/agg/$0_$1.log 2>&1 && echo "AGG OK $0 $1" || echo "AGG FAIL $0 $1"'
echo AGG DONE
