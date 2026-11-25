#!/bin/bash
# Which original predictions/optimized files contain scores >= 1.0 (DeepGO-SE CC)?
D=/ibex/scratch/projects/c2014/asaad/swissprot_proteomes_folds/deepgo-se_results/cc
for s in predictions optimized optimized_without_taxons; do
  echo "== $s"
  for f in $D/$s/*.tsv; do n=$(grep -o '|1\.0*[[:space:]]\||1\.0*$\||1[[:space:]]\||1$' $f | wc -l); [ "$n" -gt 0 ] && echo "$(basename $f) $n $(stat -c %y $f | cut -d' ' -f1)"; done
done
echo DONE
