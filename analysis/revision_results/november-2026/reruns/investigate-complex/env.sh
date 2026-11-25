source /ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/todo-reruns-20261031/common.env
IC=$BASE/todo-reruns-20261031/investigate-complex
OI=$BASE/out/investigate-complex
CPIDX=$BASE/refs/complex_validation/subunit_index.json
IDMAP=$ACCEPTED/protein-id-map.json
TAXA=(10116 1073375 1220924 187420 208964 221103 227321 229533 237561 243274 267377 272844 3197 3218 330879 331117 336722 3572 367110 3750 3847 3871 3880 3983 39947 4072 4081 510516 5270 6239 655827 69014 70601 7070 7176 7227 742740 759272 7757 7868 7955 8090 83332 8355 8364 85962 88036 9031 9913 9940)
PREDS=(seq_sim_timeset deepgo-se sprof-go)
mkdir -p $OI/logs
