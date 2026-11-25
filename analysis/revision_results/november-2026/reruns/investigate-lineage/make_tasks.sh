#!/bin/bash
# Build the termstats task list: predictor aspect taxon tau prediction-file
source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
P=$BASE/baseline-rebuild-20261029/predictions
TS=$(ls $P/seq_sim_timeset_results/cc/predictions | sed -E 's/.*taxon_([0-9]+).tsv/\1/')
HO=$(ls $P/seq_sim_heldout_results/cc/predictions | sed -E 's/.*taxon_([0-9]+).tsv/\1/')
{
for pred in seq_sim_timeset deepgo-se sprof-go; do for t in $TS; do for asp in cc mf bp; do
  tau=0.001; [ $pred = deepgo-se ] && tau=0.3
  echo $pred $asp $t $tau $(ls $P/${pred}_results/$asp/predictions/predictions_fold_*_taxon_$t.tsv)
done; done; done
for pred in seq_sim_heldout mlp_heldout; do for t in $HO; do for asp in cc mf bp; do
  echo $pred $asp $t 0.001 $(ls $P/${pred}_results/$asp/predictions/predictions_fold_*_taxon_$t.tsv)
done; done; done
} > $BASE/todo-reruns-20261031/investigate-lineage/termstats-tasks.txt
wc -l $BASE/todo-reruns-20261031/investigate-lineage/termstats-tasks.txt
