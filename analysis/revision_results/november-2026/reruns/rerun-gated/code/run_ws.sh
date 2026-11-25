#!/bin/bash
# Rerun Stage 1 + Stage 2 (standard, demotion-only, evidence-gated) for 50 timeset CC organisms on ws.
# Mirrors code-canonical-r12 slurm/revision_audit.slurm options for timeset CC, minus --sweep/--soft.
# usage: run_ws.sh NPAR [task-file]
set -uo pipefail
source /mnt/data1/leechuck/genoadjust-ws/ws.env
export R=$W/rerun-gated PFP_CODE=$CODE PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 PY CODE RES W
export TMPDIR=$R/scratch; mkdir -p $TMPDIR $R/logs $R/out
NP=${1:-12}; TASKS=${2:-$R/tasks.txt}
one() {
  p=$1; t=$2; o=$R/out/$p/cc/$t.json.gz
  [ -s $o ] && exit 0
  cd $CODE
  nice -n 10 /usr/bin/time -v $PY $R/code/revision_audit_run_gated.py --data $R/data --resources $RES \
    --predictor $p --aspect cc --taxon $t --output $o --infer-taxons --stage2 --stage2-gated \
    --gated-evidence "$W/investigate-complex/data/acc/{predictor}/cc/stage1" \
    --proteins $CODE/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv \
    --funtaxis $R/funtaxis/manifest.json --solver-seconds 120 > $R/logs/${p}_$t.log 2>&1 \
    && echo "DONE $p $t" || echo "FAIL $p $t"
}
export -f one
xargs -P $NP -n 2 bash -c 'one "$@"' _ < $TASKS
echo ALLDONE
