#!/bin/bash
# Local (lc-dell) execution of the investigation, used because IBEX per-user CPU
# quota was saturated by other jobs. Same OR-Tools (9.15.6755) and GO (sha 923645...).
# usage: run_local.sh diag|variants|floor  [P]
set -euo pipefail
export L=${L:-$HOME/tmp/claude-1000/-home-leechuck-Documents-papers-asaad-genomescale-paper/05b500ec-65d5-428c-a876-b2e360f0b256/scratchpad/icdata}
export CODE=$HOME/Documents/papers/asaad-genomescale-paper/genome-scale-pfp-adjust
export IC=$(cd $(dirname $0); pwd)
export PY=$HOME/miniconda/bin/python
export PYTHONPATH=$CODE:$CODE/pipeline PYTHONHASHSEED=0 OMP_NUM_THREADS=1
TAXA="10116 1073375 1220924 187420 208964 221103 227321 229533 237561 243274 267377 272844 3197 3218 330879 331117 336722 3572 367110 3750 3847 3871 3880 3983 39947 4072 4081 510516 5270 6239 655827 69014 70601 7070 7176 7227 742740 759272 7757 7868 7955 8090 83332 8355 8364 85962 88036 9031 9913 9940"
MODE=$1; NP=${2:-4}
mkdir -p $L/logs
tasks() { for p in seq_sim_timeset deepgo-se sprof-go; do for t in $TAXA; do echo $p $t; done; done; }
one() {
  p=$1; t=$2; cd $IC
  case $MODE in
    diag) [ -s $L/out/diag/$p/$t.json.gz ] || nice -n 19 $PY -u diag.py --code $CODE --go $L/res/go-basic.obo --accepted $L/acc --cp_index $L/res/subunit_index.json \
        --uniprot_dir $L/uniprot_index --idmap $L/res/protein-id-map.json --predictor $p --taxon $t --out $L/out/diag/$p/$t.json.gz > $L/logs/diag_${p}_$t.log 2>&1 ;;
    variants) [ -s $L/out/variants/$p/$t.json.gz ] || nice -n 19 $PY -u variants.py --code $CODE --go $L/res/go-basic.obo --accepted $L/acc --predictor $p --taxon $t \
        --out $L/out/variants --workers ${W:-3} > $L/logs/var_${p}_$t.log 2>&1 ;;
    floor) [ "$p" != deepgo-se ] || [ -s $L/out/variants_floor015/$p/$t.json.gz ] || nice -n 19 $PY -u variants.py --code $CODE --go $L/res/go-basic.obo --accepted $L/acc --predictor $p --taxon $t \
        --out $L/out/variants_floor015 --dense-floor 0.15 --variants ev_hard,ev_prefer,ev_abstain,sx_hard,sx_prefer,sx_abstain --workers ${W:-3} > $L/logs/floor_${p}_$t.log 2>&1 ;;
  esac || echo "FAIL $MODE $p $t"
}
export -f one; export MODE W
tasks | xargs -P $NP -n 2 bash -c 'one "$@"' _
echo "ALLDONE $MODE"
