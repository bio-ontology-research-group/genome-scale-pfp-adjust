# Corrected local baselines

The inherited training/transfer labels contain propagation through relations that do not imply GO annotation inheritance. Re-evaluating old predictions with corrected test references does not correct their training labels. The local Seq-Sim and MLP-ESM2 baselines must be regenerated for the corrected-training comparison. DeepGO-SE and SPROF-GO remain supplied pretrained predictor outputs; no new claim of training-data independence is made for them.

Use the completed source-reference build described in README.md. Preserve the original train/validation protein membership and organism split files. The preparation step replaces labels from source evidence, drops proteins with no remaining reference labels, canonicalizes the inherited target vocabulary and adds valid ancestors. It rejects train/validation duplication and overlap with test taxa. It does not establish a universal sequence-identity exclusion.

```bash
python data_prep/prepare_baseline_inputs.py \
  --references /path/to/rebuilt-reference-annotations \
  --go /path/to/go-basic.obo \
  --train /path/to/original/train_proteins.tsv \
  --validation /path/to/original/val_proteins.tsv \
  --test-taxa splits/heldout/test_organisms.txt \
  --terms-root /path/to/original/data \
  --output /path/to/new-baseline-inputs
```

`training-manifest.json` records the original and corrected inputs, counts and hashes. No original split or reference file is overwritten. In the audited corpus, 83,249 training and 9,254 validation proteins retain labels; 150 and 15, respectively, lose all valid reference labels. Corrected target vocabularies contain 2,668 CC, 6,575 MF and 18,762 BP terms, including roots. Roots are excluded from the downstream accuracy metrics.

Install the recorded training environment using `requirements-training.txt` (PyTorch 2.5.1+cu121); record the GPU and driver separately. On a compute node, run each GO aspect with an independent output directory:

```bash
python baselines/mlp_esm2/train_mlp_heldout.py \
  --subontology cc \
  --train-file /path/to/new-baseline-inputs/train_proteins.tsv \
  --val-file /path/to/new-baseline-inputs/val_proteins.tsv \
  --terms-file /path/to/new-baseline-inputs/cc/terms.pkl \
  --test-organisms-file splits/heldout/test_organisms.txt \
  --go-obo /path/to/go-basic.obo \
  --esm2-embeddings-dir /path/to/embedding-parent \
  --output-dir /path/to/new-predictions/mlp_heldout_results/cc \
  --device cuda:0 --epochs 10 --batch-size 256 --lr 0.001 \
  --seed 42 --sparsity-threshold 0.001 --num-workers 4
```

The embedding parent contains `esm2_embeddings/<taxon>/esm2.pkl`. Missing, duplicate, nonfinite or incorrectly sized features are errors. No zero vectors are fabricated and dimensions are not silently padded or truncated. Training uses CPU-resident data with GPU minibatches; inference validates every test organism. The model checkpoint is saved before validation/inference. Result metadata include embedding hashes, software, training parameters and checkpoint identity. Validation metrics in the training log are micro-averaged diagnostics, not the protein-macro metrics used in the paper. The legacy `--no-save-predictions` / `--annotations-go-basic-dir` evaluation path is retired: export all test predictions, then evaluate them with the shared revision pipeline and source-rebuilt reference labels.

Seq-Sim rebuilds both sequence databases and scores from the verified source DAT files and corrected reference labels. Run once per organism split; all three aspects are generated together:

```bash
python baselines/seq_sim/rebuild_from_sources.py \
  --references /path/to/rebuilt-reference-annotations \
  --source-root /path/to/uniprot_reference_proteomes \
  --go /path/to/go-basic.obo \
  --terms-root /path/to/new-baseline-inputs \
  --train-taxa /path/to/train_organisms.txt \
  --test-taxa /path/to/test_organisms.txt \
  --output /path/to/new-predictions/seq_sim_timeset_results \
  --diamond /path/to/diamond --threads 16 --score-floor 0.001
```

DIAMOND uses sensitive mode, E-value 1e-5, at most 50 target sequences and minimum bit score 50. Transfer scores are the fraction of total matched-subject bit score assigned to a term. Every queried source protein is written, including empty rows for proteins without hits. The prediction manifest records input hashes, executable version, exact commands, complete query counts and output hashes. The old `predict_seq_similarity.py` CLI is retired because it accepted unverified propagated labels and omitted proteins without aspect-specific predictions.

Both local baselines serialize scores at 17-digit precision and omit scores below 0.001. Downstream repair at the strict threshold 0.001 therefore preserves all active pairs; continuous-score accuracy is explicitly the accuracy of these serialized scores, with omitted scores zero. This is not a claim to evaluate untruncated dense model outputs. Apply the same canonical alignment and maximum-ancestor propagation as the other predictors before comparisons.

```bash
python tests/test_baseline_inputs.py
python tests/test_sequence_transfer.py
python tests/test_mlp_contracts.py  # needs the training/PyTorch environment
```
