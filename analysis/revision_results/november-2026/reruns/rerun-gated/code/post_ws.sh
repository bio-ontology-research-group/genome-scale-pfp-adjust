#!/bin/bash
# Post-processing on ws: reproduction check, aggregation (paper evaluator), Complex Portal validation, UniProt support.
set -euo pipefail
source /mnt/data1/leechuck/genoadjust-ws/ws.env
R=$W/rerun-gated; L=$W/investigate-complex/data; IC=$W/investigate-complex
export PFP_CODE=$CODE TMPDIR=$R/scratch OMP_NUM_THREADS=1
mkdir -p $R/summaries $R/membership $R/eval/vdir
cd $R
:
for p in seq_sim_timeset deepgo-se sprof-go; do
  nice -n 10 $PY code/aggregate_revision_audit_gated.py --tasks tasks3.txt --results out --go $RES/go-basic.obo --annotations $REFS \
    --proteins $CODE/splits/timeset/proteins_by_date_23-MAY-2024_filtered.tsv --heldout-train $CODE/splits/heldout/train_organisms.txt \
    --predictor $p --aspect cc --output summaries/${p}_cc.json > logs/aggregate_$p.log 2>&1 &
done
for p in seq_sim_timeset deepgo-se sprof-go; do tau=0.001; [ $p = deepgo-se ] && tau=0.3
  for m in; do
    nice -n 10 $PY $CODE/analysis/complex_validation/validate_stage2.py --stage1_dir out/$p/cc/stage1 --stage2_dir out/$p/cc/stage2_$m \
      --subunit_index $L/res/subunit_index.json --protein-id-map $L/res/protein-id-map.json --go $RES/go-basic.obo \
      --threshold $tau --label "$p $m" --output_file membership/membership_${p}_$m.json > logs/membership_${p}_$m.log 2>&1 &
  done
  mkdir -p eval/vdir/$p; for t in $(cat taxa.txt); do ln -sf $R/out/$p/cc/gated_variants/$t.json.gz eval/vdir/$p/$t.json.gz; done
done
wait
cd $IC/scripts && nice -n 10 $PY eval_variants.py --code $CODE --go $RES/go-basic.obo --accepted $R/out --cp_index $L/res/subunit_index.json \
  --uniprot_dir $L/uniprot_index --idmap $L/res/protein-id-map.json --vdir $R/eval/vdir --taxa $(paste -sd, $R/taxa.txt) --out $R/eval/eval_gated.json > $R/logs/eval.log 2>&1
echo POSTDONE
