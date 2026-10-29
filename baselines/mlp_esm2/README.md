# MLP-ESM2 baseline

The corrected baseline trains an MLP for each GO aspect on supplied 2,560-dimensional ESM-2 embeddings and predicts all proteins in the 219 heldout organisms. Its hidden and residual blocks have 1,024 units, with ReLU, layer normalization and dropout 0.1; a sigmoid output produces term scores.

Follow [the corrected baseline workflow](../../reproducibility/BASELINES.md) for source-reference rebuilding, verified train/validation inputs, the pinned PyTorch environment and complete commands. The historical SLURM training wrapper is retired. Original propagated annotation files must not be reused as corrected training labels.

`train_mlp_heldout.py` requires the training-input manifest, rejects missing or malformed features, records consumed embedding hashes and saves the fitted checkpoint before validation and prediction export. Validation metrics are micro-averaged training diagnostics. Paper test metrics are computed from serialized predictions by the shared revision evaluator; the legacy annotated-only test-evaluation option is retired.

`extract_esm2_embeddings.py` is the historical extraction utility for mean-pooled `esm2_t36_3B_UR50D` features. The corrected experiment uses the supplied, checksum-recorded embeddings; it does not claim to have regenerated the pretrained representation. The training CLI's embedding parent contains `esm2_embeddings/<taxon>/esm2.pkl`.

Reference: Lin et al. (2023), *Science* 379:1123–1130 (ESM-2 backbone).
