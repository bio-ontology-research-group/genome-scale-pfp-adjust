#!/bin/bash
R=/mnt/data1/leechuck/genoadjust-ws/rerun-notaxon; cd $R
scripts/run_tasks.sh tasks_heldout.txt 1 > logs/run_heldout_retry2.log 2>&1
scripts/aggregate.sh tasks_heldout.txt 1 mlp_heldout > logs/agg_mlp.log 2>&1
echo CHAIN3 DONE
