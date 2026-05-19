# SLURM job scripts

All scripts source `config.env` at the top, so values like `DATA_DIR`,
`CONDA_BASE`, `CONDA_ENV`, `GAEF_DIR`, `GO_OBO` resolve consistently.

## Quickstart

```
cp config.env.example config.env
$EDITOR config.env                     # fill in paths
sbatch slurm/adjust/seq_sim_timeset.slurm
```

## Layout

| Phase | Script | Args | What it does |
|---|---|---|---|
| Embed   | `train/extract_esm2.slurm` | — | ESM2 mean embeddings per proteome (needed by MLP-ESM2 + DeepGO-SE) |
| Predict | `train/predict_seq_sim_heldout.slurm` | — | Diamond-BLAST Seq-Sim baseline on the heldout split |
| Predict | `train/predict_seq_sim_timeset.slurm` | — | Diamond-BLAST Seq-Sim baseline on the timeset split |
| Predict | `train/train_mlp_heldout.slurm`       | `TASK` (`cc`/`mf`/`bp`) | Train MLP-ESM2 on the heldout train/val partition and predict on its test organisms |
| Predict | `train/predict_deepgo_se.slurm`       | `TASK` (`cc`/`mf`/`bp`) | Pretrained DeepGO-SE ensemble inference on the timeset test set |
| Predict | `train/predict_sprof_go.slurm`        | — | Pretrained SPROF-GO inference on the timeset test set |
| Adjust  | `adjust/seq_sim_timeset.slurm` | — | Seq-Sim (timeset): Stage 1 with/without taxon + Stage 2 on CC |
| Adjust  | `adjust/seq_sim_heldout.slurm` | — | Seq-Sim (heldout): Stage 1 with/without taxon |
| Adjust  | `adjust/mlp_heldout.slurm`     | `[both\|with-taxons\|without-taxons]` | MLP-ESM2 (heldout): Stage 1 |
| Adjust  | `adjust/deepgo_se.slurm`       | — | DeepGO-SE (timeset): Stage 1 with/without taxon + Stage 2 on CC |
| Adjust  | `adjust/sprof_go.slurm`        | `TASK` (`cc`/`mf`/`bp`) | SPROF-GO (timeset): Stage 1 with/without taxon, plus Stage 2 when TASK=cc |
| Eval    | `eval/gaef_annotations.slurm` | — | GAEF metrics on ground-truth annotations |
| Eval    | `eval/gaef_timeset_methods.slurm` | — | GAEF metrics on each predictor (array, 9 tasks) |
| Ablation | `ablation/heuristic.slurm` | `TAXON` | H1/H2/H3 ablation on one test organism |

## Required env vars

Set in `config.env` (or override on the `sbatch` line):

- `DATA_DIR` — root of working data (predictions, annotations, intermediate)
- `CONDA_BASE`, `CONDA_ENV` — conda setup (env is built from `environment.yml`)
- `GAEF_DIR` — path to the GAEF sibling repo
- `DEEPGO_SE_REPO` — path to the deepgo2 clone (read by `train/predict_deepgo_se.slurm`)
- `SPROF_GO_REPO` — path to the SPROF-GO clone (read by `train/predict_sprof_go.slurm`)
- `GO_OBO` — path to `go-basic.obo` (defaults to `data/go-basic.obo`)

Optional per-script overrides shown in each script's header comment.
