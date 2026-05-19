# Dataset splits

Two splits used in the paper. Organism lists are committed; per-protein TSVs
are regenerated from the organism lists + Swiss-Prot via `data_prep/`.

## `timeset/` — 50 organisms

Test proteins are those with experimental Swiss-Prot integration date strictly
after **2024-05-23**, chosen to fall after the training cutoff of the pretrained
predictors (DeepGO-SE, SPROF-GO). The rest of the proteomes form the training
pool.

- `test_organisms.txt` — 50 taxa containing at least one post-cutoff protein
- `train_organisms.txt` — remaining taxa
- `proteins_by_date_23-MAY-2024_filtered.tsv` — the cutoff-filtered protein
  table; the source for the timeset test proteins

Regenerate via `python data_prep/extract_proteins_by_date.py --cutoff-date 23-MAY-2024 ...`.

## `heldout/` — 219 organisms

Test set is sampled to contain ~10% of all experimental annotations across the
benchmark. Train/val are partitioned 90:10 using MMseqs2 easy-cluster
(`--min-seq-id 0.3`): no two proteins with ≥30% sequence identity cross the
train/val boundary.

- `test_organisms.txt` — 219 held-out taxa
- `train_organisms.txt` — taxa contributing to the train/val pool
- `split_info.json` — split summary (counts, seed, parameters)

Regenerate via `python data_prep/create_heldout_split.py ...` after the
benchmark TSVs have been built by `data_prep/build_dataset.py`.
