#!/bin/bash
# Usage: run_tasks.sh TASKFILE NPAR   (TASKFILE lines: predictor aspect taxon)
source /mnt/data1/leechuck/genoadjust-ws/ws.env
export PYTHONPATH=$PYTHONPATH:$W/rerun-notaxon; export CODE PYTHONPATH PYTHONHASHSEED
R=$W/rerun-notaxon; export R RES REFS PY
TASKS=$1; NPAR=$2
while read -r pred asp t; do
  [ -s $R/out/$pred/$asp/$t.json.gz ] || echo "$pred $asp $t"
done < $TASKS | xargs -P $NPAR -L 1 sh -c '
  pred=$0; asp=$1; t=$2
  case $pred in *_heldout) lvd=$R/../investigate-lineage/lv-heldout; extra="";;
     *) lvd=$R/../investigate-lineage/lv; extra="--proteins $CODE/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv";; esac
  mkdir -p $R/out/$pred/$asp $R/logs/run
  s=$(date +%s)
  if nice -n 10 /usr/bin/time -v $PY -u $R/scripts/notaxon_run.py --data $R/data --resources $RES --annotations $REFS \
     --predictor $pred --aspect $asp --taxon $t --output $R/out/$pred/$asp/$t.json.gz --infer-taxons \
     --lineage-dir $lvd $extra > $R/logs/run/${pred}_${asp}_$t.log 2>&1; then
    echo "OK $pred $asp $t $(( $(date +%s)-s ))s $(grep Maximum.resident $R/logs/run/${pred}_${asp}_$t.log | awk "{print \$NF}")KB"
  else echo "FAIL $pred $asp $t"; fi'
echo ALL DONE $TASKS
