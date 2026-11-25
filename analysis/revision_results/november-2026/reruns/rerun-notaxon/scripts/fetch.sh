#!/bin/bash
# fetch prediction files for predictors matching regex $1 into $W/data/<pred>/<asp>/ (the layout <pred>_results/<asp>/predictions/)
source /mnt/data1/leechuck/genoadjust-ws/ws.env
R=$W/rerun-notaxon; D=$R/data
awk -v re="$1" '$1 ~ re {print $1, $2, $5}' $W/investigate-lineage/termstats-tasks.txt | while read -r pred asp f; do
  d=$D/${pred}_results/$asp/predictions; mkdir -p $d
  [ -s $d/$(basename $f) ] || echo "$f $d"
done | xargs -P 4 -L 1 sh -c 'rsync -a ibex:$0 $1/ 2>/dev/null || echo FAIL $0'
echo FETCH DONE $1
