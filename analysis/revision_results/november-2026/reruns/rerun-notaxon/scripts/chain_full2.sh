#!/bin/bash
R=/mnt/data1/leechuck/genoadjust-ws/rerun-notaxon; cd $R
NEED_GB=24 scripts/aggregate_full.sh tasks_timeset.txt 3 seq_sim_timeset deepgo-se sprof-go > logs/agg_full_timeset2.log 2>&1
grep -v "mlp_heldout bp" tasks_heldout.txt > tasks_heldout_nomlpbp.txt; grep "mlp_heldout bp" tasks_heldout.txt > tasks_heldout_mlpbp.txt
NEED_GB=16 scripts/run_tasks_full.sh tasks_heldout_nomlpbp.txt 6 > logs/run_full_heldout2.log 2>&1
NEED_GB=40 scripts/run_tasks_full.sh tasks_heldout_mlpbp.txt 3 > logs/run_full_heldout_mlpbp2.log 2>&1
NEED_GB=40 scripts/run_tasks_full.sh tasks_heldout.txt 2 > logs/run_full_heldout_retry2.log 2>&1
NEED_GB=48 scripts/aggregate_full.sh tasks_heldout.txt 2 seq_sim_heldout mlp_heldout > logs/agg_full_heldout2.log 2>&1
echo CHAIN FULL2 DONE
