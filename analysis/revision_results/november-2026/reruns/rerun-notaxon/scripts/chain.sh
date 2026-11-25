#!/bin/bash
R=/mnt/data1/leechuck/genoadjust-ws/rerun-notaxon; cd $R
until grep -q "FETCH DONE" logs/fetch_timeset.log; do sleep 30; done
scripts/run_tasks.sh tasks_timeset.txt 16 > logs/run_timeset.log 2>&1
until grep -q "FETCH DONE" logs/fetch_heldout.log; do sleep 30; done
scripts/run_tasks.sh tasks_heldout.txt 12 > logs/run_heldout.log 2>&1
