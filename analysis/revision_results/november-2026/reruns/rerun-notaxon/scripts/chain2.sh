#!/bin/bash
R=/mnt/data1/leechuck/genoadjust-ws/rerun-notaxon; cd $R
scripts/run_tasks.sh tasks_timeset.txt 4 > logs/run_timeset_retry.log 2>&1
scripts/run_tasks.sh tasks_heldout.txt 8 > logs/run_heldout2.log 2>&1
scripts/run_tasks.sh tasks_heldout.txt 3 > logs/run_heldout_retry.log 2>&1
echo CHAIN2 DONE
