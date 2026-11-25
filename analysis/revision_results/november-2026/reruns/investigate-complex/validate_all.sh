#!/bin/bash
# Run the unchanged Complex Portal validator on every variant, then the summary.
set -euo pipefail
L=${L:-$HOME/tmp/claude-1000/-home-leechuck-Documents-papers-asaad-genomescale-paper/05b500ec-65d5-428c-a876-b2e360f0b256/scratchpad/icdata}
CODE=$HOME/Documents/papers/asaad-genomescale-paper/genome-scale-pfp-adjust
IC=$(cd $(dirname $0); pwd); PY=$HOME/miniconda/bin/python
export PYTHONPATH=$CODE:$CODE/pipeline PYTHONHASHSEED=0
mkdir -p $L/out/membership
for vd in $L/out/variants $L/out/variants_floor015; do
  for pdir in $vd/*/; do p=$(basename $pdir); tau=0.001; [ $p = deepgo-se ] && tau=0.3
    for v in $pdir*/; do v=$(basename $v); [ $(ls $pdir$v | wc -l) -eq 50 ] || { echo "incomplete $p $v"; continue; }
      nice -n 19 $PY $CODE/analysis/complex_validation/validate_stage2.py --stage1_dir $L/acc/$p/cc/stage1 --stage2_dir $pdir$v \
        --subunit_index $L/res/subunit_index.json --protein-id-map $L/res/protein-id-map.json --go $L/res/go-basic.obo \
        --threshold $tau --label "$p $v" --output_file $L/out/membership/$(basename $vd)_${p}_$v.json > /dev/null &
    done; wait
  done
done
TAXA=$(ls $L/acc/deepgo-se/cc/stage1 | sed 's/.*taxon_//;s/.tsv//' | paste -sd,)
cd $IC
for vd in variants variants_floor015; do
  nice -n 19 $PY eval_variants.py --code $CODE --go $L/res/go-basic.obo --accepted $L/acc --cp_index $L/res/subunit_index.json \
    --uniprot_dir $L/uniprot_index --idmap $L/res/protein-id-map.json --vdir $L/out/$vd --taxa $TAXA --out $L/out/eval_$vd.json > /dev/null
done
echo VALIDATED
