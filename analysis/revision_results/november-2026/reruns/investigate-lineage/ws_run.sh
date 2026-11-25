#!/bin/bash
# ws driver (no Slurm). Usage:
#   ws_run.sh termstats NPAR            # fetch + summarize termstats missing locally
#   ws_run.sh variants VARS SUB REGEX NPAR   # lineage_variants for units of predictors matching REGEX
#   ws_run.sh contigs VARS SUB NPAR     # contig experiment, 15 combos
set -u
source /mnt/data1/leechuck/genoadjust-ws/ws.env
IL=$W/investigate-lineage; S=$IL/scripts; TS=$IL/termstats; mkdir -p $IL/logs
mode=$1; shift
case $mode in
termstats)
  NPAR=$1
  while read -r pred asp t tau f; do
    [ -s $TS/$pred/$asp/$t.json.gz ] || echo "$pred $asp $t $tau $f"
  done < $IL/termstats-tasks.txt | xargs -P $NPAR -L 1 sh -c '
    pred=$0; asp=$1; t=$2; tau=$3; f=$4; d='$W'/data/$pred/$asp; mkdir -p $d '$TS'/$pred/$asp
    rsync -a ibex:$f $d/ 2>/dev/null
    nice -n 10 '$PY' '$S'/termstats.py --code '$CODE' --resources '$RES' --predictions $d/$(basename $f) --aspect $asp --tau $tau --output '$TS'/$pred/$asp/$t.json.gz > /dev/null && echo "$pred $asp $t ok"
    rm -f $d/$(basename $f)'
  ;;
variants)
  VARS=$1; SUB=$2; RE=$3; NPAR=$4; mkdir -p $IL/$SUB
  awk -v re="$RE" '$2=="cc" && $1 ~ re {print $1, $3}' $IL/termstats-tasks.txt | while read -r pred t; do
    [ -e $IL/$SUB/.done_${pred}_$t ] && continue
    [ -s $TS/$pred/cc/$t.json.gz ] && [ -s $TS/$pred/mf/$t.json.gz ] && [ -s $TS/$pred/bp/$t.json.gz ] && echo "$pred $t"
  done | xargs -P $NPAR -L 1 sh -c '
    nice -n 10 '$PY' -u '$S'/lineage_variants.py --code '$CODE' --resources '$RES' --termstats '$TS' --predictor $0 --taxa $1 \
      --variants '$VARS' --truth --ic-file '$IL'/go_ic.tsv --time-limit 120 --outdir '$IL/$SUB' 2>&1 | grep -v -E "^\[INFO\]|^Computed" > '$IL'/logs/'$SUB'_$0_$1.log \
      && touch '$IL/$SUB'/.done_$0_$1'
  ;;
contigs)
  VARS=$1; SUB=$2; NPAR=$3; K=${4:-1000}; mkdir -p $IL/$SUB $W/data/contigs
  D=$W/todo-reruns/contigs
  for T in 83332 3572 69014; do for c in "seq_sim_timeset bp" "deepgo-se bp" "sprof-go bp" "seq_sim_timeset cc" "seq_sim_timeset mf"; do echo $T $c; done; done |
  xargs -P $NPAR -L 1 sh -c '
    T=$0; PRED=$1; ASP=$2; TAU=0.001; [ $PRED = deepgo-se ] && TAU=0.3
    f=$(awk -v p=$PRED -v a=$ASP -v t=$T "\$1==p && \$2==a && \$3==t {print \$5}" '$IL'/termstats-tasks.txt)
    loc='$W'/data/contigs/${T}_${PRED}_${ASP}.tsv; [ -s $loc ] || rsync -a ibex:$f $loc
    SIZE=$(grep "Total concatenated" '$D'/cds_$T.log | sed -E "s/.*size: ([0-9,]+) bp.*/\1/" | tr -d ,)
    nice -n 10 '$PY' -u '$S'/contig_variants.py --code '$CODE' --resources '$RES' --predictions $loc --cds '$D'/cds_$T.tsv \
      --genome-size $SIZE --taxon $T --aspect $ASP --tau $TAU --variants '$VARS' --evaluate-first '$K' --output '$IL/$SUB'/${T}_${PRED}_${ASP}.json \
      2>&1 | grep -v -E "^\[INFO\]|^Computed" > '$IL'/logs/'$SUB'_${T}_${PRED}_${ASP}.log'
  ;;
esac
