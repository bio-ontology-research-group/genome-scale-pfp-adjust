"""
Script to train MLP on the UniProt proteomes using a fixed train/val/test split.
Save predictions for each test organism in the format: protein_id<tab>GO:term|score<tab>...

Given:
- train_proteins.tsv: protein_id<tab>taxon_id<tab>GO:term1<tab>GO:term2<tab>...
- val_proteins.tsv:   protein_id<tab>taxon_id<tab>GO:term1<tab>GO:term2<tab>...
- test_organisms.txt: one taxon_id per line
- esm2_embeddings_dir directory:
    - contains the esm2.pkl file for each taxon
- output_dir directory:
    - predictions/ subdir with predictions_fold_XX_taxon_YYYY.tsv (protein_id<tab>GO:term|score<tab>...)

"""

import click as ck
import os
import json
from typing import Dict, List, Tuple, Set
from collections import defaultdict
import time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
import gc
import shutil
import platform
from multiprocessing import Pool

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))  # to have access to deepgo modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'GAEF')))  # to have access to GAEF modules

# NumPy compatibility shim: pickle files produced by NumPy 2.x use numpy._core,
# but older NumPy versions expose the same symbols under numpy.core.
try:
    import numpy._core.numeric  # type: ignore # noqa: F401
except ModuleNotFoundError:
    import numpy.core as _np_core
    sys.modules['numpy._core'] = _np_core
    sys.modules['numpy._core.numeric'] = _np_core.numeric
    sys.modules['numpy._core.multiarray'] = _np_core.multiarray
    sys.modules['numpy._core._multiarray_umath'] = getattr(_np_core, '_multiarray_umath', _np_core.multiarray)

from sklearn.metrics import roc_auc_score
from deepgo.base import Residual, MLPBlock
from gaef_patches.ontology import Ontology
from data_prep.prepare_baseline_inputs import digest

EMBEDDING_INPUTS = {}

# Route click echo through tqdm to avoid breaking progress bars
try:
    _original_ck_echo = ck.echo
    def _tqdm_echo(message: str = "", file=None, nl: bool = True, err: bool = False, color=None):
        tqdm.write(str(message))
    ck.echo = _tqdm_echo  # type: ignore
except Exception:
    pass


import random


def set_seed(seed: int = 42):
    """Seed all RNGs (Python, NumPy, PyTorch CPU/CUDA) for reproducible training.

    Added 2026-06 so the reported MLP-ESM2 numbers can be reproduced from a fixed
    seed; CP-SAT downstream is already deterministic given its input. See the paper
    Methods (Baseline predictors).
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


class MLPModel(nn.Module):
    """
    Baseline MLP model with two fully connected layers with residual connection.
    """

    def __init__(self, input_length, nb_gos, nodes=[1024,]):
        super().__init__()
        self.nb_gos = nb_gos
        net = []
        for hidden_dim in nodes:
            net.append(MLPBlock(input_length, hidden_dim))
            net.append(Residual(MLPBlock(hidden_dim, hidden_dim)))
            input_length = hidden_dim
        net.append(nn.Linear(input_length, nb_gos))
        net.append(nn.Sigmoid())
        self.net = nn.Sequential(*net)

    def forward(self, features):
        features = features.reshape(features.shape[0], -1)
        return self.net(features)


def load_proteins_from_split_tsv(filepath: str) -> pd.DataFrame:
    """
    Load proteins from a train/val split TSV file.

    File format per line: protein_id<tab>taxon_id<tab>GO:term1<tab>GO:term2<tab>...

    Args:
        filepath: Path to the TSV file (train_proteins.tsv or val_proteins.tsv)

    Returns:
        DataFrame with columns: proteins, prop_annotations, orgs
    """
    all_proteins = []
    all_annotations = []
    all_orgs = []

    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            protein_id = parts[0].strip()
            taxon_id = parts[1].strip()
            go_terms = [t.strip() for t in parts[2:] if t.strip()]
            if protein_id and taxon_id:
                all_proteins.append(protein_id)
                all_annotations.append(go_terms)
                all_orgs.append(taxon_id)

    df = pd.DataFrame({
        'proteins': all_proteins,
        'prop_annotations': all_annotations,
        'orgs': all_orgs,
    })
    return df


def load_esm2_embeddings_taxon(esm2_embeddings_dir: str, taxon_id: str) -> Dict[str, np.ndarray]:
    """
    Load ESM2 embeddings from pickle file for a given taxon.

    Args:
        esm2_embeddings_dir: Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl
        taxon_id: Taxon ID

    Returns:
        Dictionary mapping protein_id to numpy array (embedding vector)
    """
    esm2_file = os.path.join(esm2_embeddings_dir, 'esm2_embeddings', taxon_id, 'esm2.pkl')
    if not os.path.exists(esm2_file):
        raise FileNotFoundError(esm2_file)
    EMBEDDING_INPUTS[os.path.abspath(esm2_file)] = digest(esm2_file)
    df = pd.read_pickle(esm2_file)
    if not {'protein', 'esm2'}.issubset(df.columns):
        raise ValueError(f"Missing protein/esm2 columns: {esm2_file}")
    embeddings = {}
    for row in df.itertuples():
        raw_id = str(row.protein)
        protein_id = raw_id.split('|')[2].split()[0] if '|' in raw_id else raw_id.split()[0]
        if protein_id in embeddings:
            raise ValueError(f"Duplicate embedding ID: {protein_id}")
        vec = np.asarray(row.esm2, dtype=np.float32)
        if vec.ndim != 1 or not np.isfinite(vec).all():
            raise ValueError(f"Invalid embedding: {protein_id}")
        embeddings[protein_id] = vec
    if not embeddings:
        raise ValueError(f"Empty embedding file: {esm2_file}")
    return embeddings


def load_features_uniprot(df: pd.DataFrame, esm2_embeddings_dir: str, vec_dim: int) -> Tuple[torch.Tensor, List[str]]:
    """
    Load ESM2 embeddings for proteins in DataFrame.

    Args:
        df: DataFrame with 'proteins' and 'orgs' columns
        esm2_embeddings_dir: Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl files
        vec_dim: Expected vector dimension

    Returns:
        Tuple of (features_tensor [N, D], missing_proteins list)
    """
    features_np = np.zeros((len(df), vec_dim), dtype=np.float32)
    missing: List[str] = []

    # Load one proteome at a time; retaining every proteome exhausts RAM.
    rows_by_taxon = defaultdict(list)
    for i, row in enumerate(df.itertuples()):
        rows_by_taxon[str(row.orgs)].append((i, str(row.proteins)))
    for taxon_id, rows in rows_by_taxon.items():
        embeddings = load_esm2_embeddings_taxon(esm2_embeddings_dir, taxon_id)
        for i, protein_id in rows:
            if protein_id not in embeddings:
                missing.append(protein_id)
                continue
            vec = embeddings[protein_id]
            if vec.shape != (vec_dim,):
                raise ValueError(f"Embedding dimension differs for {protein_id}: {vec.shape}, expected {vec_dim}")
            features_np[i] = vec
        del embeddings
    if missing:
        raise ValueError(f"Missing embeddings for {len(missing)} proteins: {missing[:10]}")

    features = torch.from_numpy(features_np)
    return features, missing


def _prepare_labels(df: pd.DataFrame, terms: List[str]) -> torch.Tensor:
    """
    Create multi-hot label tensor [N, T] given df with 'prop_annotations'.
    """
    term_index = {t: i for i, t in enumerate(terms)}
    labels = np.zeros((len(df), len(terms)), dtype=np.float32)
    for i, row in enumerate(df.itertuples()):
        annots = row.prop_annotations if 'prop_annotations' in df.columns else []
        if annots is None or (isinstance(annots, float) and pd.isna(annots)):
            continue
        for go_id in annots:
            j = term_index.get(go_id)
            if j is not None:
                labels[i, j] = 1.0
    return torch.from_numpy(labels)


def train_mlp(train_df: pd.DataFrame, terms: List[str], device: str, esm2_embeddings_dir: str, vec_dim: int,
              epochs: int = 10, batch_size: int = 256, lr: float = 1e-3) -> Tuple[MLPModel, Dict[str, int]]:
    """
    Train MLP on provided training dataframe and list of GO terms.
    Returns trained model and a small info dict with counts.
    """
    X, missing = load_features_uniprot(train_df, esm2_embeddings_dir, vec_dim)
    y = _prepare_labels(train_df, terms)

    input_dim = int(X.shape[1])
    output_dim = int(len(terms))

    model = MLPModel(input_dim, output_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.BCELoss()

    ds = TensorDataset(X, y)
    dl = DataLoader(ds, batch_size=batch_size, shuffle=True, drop_last=False)

    model.train()
    for epoch in tqdm(range(epochs), desc="Training MLP", leave=True):
        epoch_loss = 0.0
        for xb, yb in dl:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            probs = model(xb)
            loss = criterion(probs, yb)
            loss.backward()
            optimizer.step()
            epoch_loss += float(loss.detach().cpu().item())
        ck.echo(f"Epoch {epoch+1}/{epochs} - loss {epoch_loss / max(1, len(dl)):.4f}")

    info = {"n_samples": int(len(train_df)), "n_missing_features": int(len(missing))}
    return model, info


def compute_fmax(labels: np.ndarray, preds_proba: np.ndarray) -> Tuple[float, float, float, float]:
    """
    Compute Fmax metric by testing different thresholds and finding the one that maximizes F1 score.
    Optimized vectorized version.

    Args:
        labels: Binary labels [N, T]
        preds_proba: Prediction probabilities [N, T]

    Returns:
        Tuple of (fmax, optimal_threshold, precision_at_fmax, recall_at_fmax)
    """
    thresholds = np.arange(0.0, 1.01, 0.01)

    labels_flat = labels.astype(bool).flatten()
    preds_flat = preds_proba.flatten()

    # Evaluate one threshold at a time instead of allocating 101 dense matrices.
    tp, fp, fn = [], [], []
    for threshold in thresholds:
        selected = (preds_flat > 0) & (preds_flat >= threshold)
        tp.append(np.count_nonzero(selected & labels_flat))
        fp.append(np.count_nonzero(selected & ~labels_flat))
        fn.append(np.count_nonzero(~selected & labels_flat))
    tp, fp, fn = np.asarray(tp), np.asarray(fp), np.asarray(fn)

    with np.errstate(divide='ignore', invalid='ignore'):
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * precision * recall / (precision + recall)

        precision = np.nan_to_num(precision, nan=0.0)
        recall = np.nan_to_num(recall, nan=0.0)
        f1 = np.nan_to_num(f1, nan=0.0)

    best_idx = np.argmax(f1)

    return float(f1[best_idx]), float(thresholds[best_idx]), \
           float(precision[best_idx]), float(recall[best_idx])


def find_optimal_threshold_on_validation(model: MLPModel, validation_df: pd.DataFrame, terms: List[str],
                                         device: str, esm2_embeddings_dir: str,
                                         vec_dim: int) -> Dict[str, float]:
    """
    Find optimal threshold using validation set.

    Returns:
        Dict with keys: optimal_threshold, fmax, precision, recall, auc
    """
    model.eval()
    with torch.no_grad():
        X_val, _ = load_features_uniprot(validation_df, esm2_embeddings_dir, vec_dim)
        y_val = _prepare_labels(validation_df, terms)
        preds_proba = np.concatenate([
            model(X_val[i:i+256].to(device)).cpu().numpy()
            for i in range(0, len(X_val), 256)
        ])
        labels = y_val.numpy()

        ck.echo(f"[INFO] Validation set: {len(validation_df)} proteins, {len(terms)} terms")
        fmax, optimal_threshold, precision, recall = compute_fmax(labels, preds_proba)

        # Micro-averaged ROC AUC (flatten to treat each label independently)
        try:
            auc = float(roc_auc_score(labels.flatten(), preds_proba.flatten()))
        except ValueError:
            auc = float('nan')

        ck.echo(f"[INFO] Validation metrics:")
        ck.echo(f"  Fmax:               {fmax:.4f}")
        ck.echo(f"  Optimal threshold:  {optimal_threshold:.4f}")
        ck.echo(f"  Precision at Fmax:  {precision:.4f}")
        ck.echo(f"  Recall at Fmax:     {recall:.4f}")
        ck.echo(f"  AUC (micro):        {auc:.4f}")

        return {
            'optimal_threshold': optimal_threshold,
            'fmax': fmax,
            'precision': precision,
            'recall': recall,
            'auc': auc,
        }


def evaluate_on_annotated_test_proteins(
    model: MLPModel,
    test_taxon_ids: List[str],
    annotations_go_basic_dir: str,
    terms: List[str],
    device: str,
    esm2_embeddings_dir: str,
    vec_dim: int,
    output_dir: str,
    batch_size: int = 1024,
) -> Dict:
    """Retired: evaluate serialized outputs through the common revision evaluator."""
    raise RuntimeError(
        "The legacy annotated-test evaluator is retired. Export predictions and use "
        "analysis/aggregate_revision_audit.py with source-rebuilt reference labels."
    )


def _propagate_single_protein(protein_id: str, go_scores: Dict[str, float],
                               ancestor_cache: Dict[str, Set[str]]) -> Dict[str, float]:
    """
    Propagate scores for a single protein through ontology hierarchy.
    """
    propagated_scores: Dict[str, float] = {}

    for go_id, score in go_scores.items():
        if go_id not in propagated_scores or propagated_scores[go_id] < score:
            propagated_scores[go_id] = score

        ancestors = ancestor_cache.get(go_id, set())
        for anc in ancestors:
            if anc not in propagated_scores or propagated_scores[anc] < score:
                propagated_scores[anc] = score

    return propagated_scores


def _process_chunk_and_write(args: Tuple) -> str:
    """
    Worker function to propagate scores and write chunk to temp file.

    Args:
        args: Tuple of (chunk_id, protein_scores_list, ancestor_cache, sparsity_threshold, temp_file)
              where protein_scores_list is List[Tuple[protein_id, Dict[go_term, score]]]

    Returns:
        Path to the temp file written
    """
    chunk_id, protein_scores_list, ancestor_cache, sparsity_threshold, temp_file = args

    with open(temp_file, 'w') as f:
        for protein_id, go_scores in protein_scores_list:
            propagated = _propagate_single_protein(protein_id, go_scores, ancestor_cache)

            if '|' in protein_id:
                protein_id = protein_id.split('|')[2].split()[0]

            if not propagated:
                f.write(f"{protein_id}\n")
            else:
                filtered_terms = [f"{go_term}|{score:.17g}"
                                  for go_term, score in propagated.items()
                                  if score >= sparsity_threshold]
                if filtered_terms:
                    f.write(f"{protein_id}\t" + "\t".join(filtered_terms) + "\n")
                else:
                    f.write(f"{protein_id}\n")

    return temp_file


def generate_propagate_write_parallel(model: nn.Module, esm2_embeddings_dir: str, taxon_id: str,
                                       terms: List[str], device: str, vec_dim: int,
                                       go: Ontology, ancestor_cache: Dict[str, Set[str]],
                                       output_file: str, sparsity_threshold: float = 0.0,
                                       batch_size: int = 1024, num_workers: int = 4,
                                       chunk_size: int = 5000) -> int:
    """
    Generate predictions, propagate scores, and write to TSV using parallel workers.

    Main process handles MLP inference (GPU-bound).
    Worker processes handle propagation and writing to temp files (CPU-bound).
    Temp files are merged at the end.

    Returns:
        Number of proteins processed
    """
    embeddings = load_esm2_embeddings_taxon(esm2_embeddings_dir, taxon_id)
    if not embeddings:
        return 0

    protein_ids = list(embeddings.keys())
    num_proteins = len(protein_ids)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Pre-fetch ancestors for all terms
    missing_terms = set(terms) - ancestor_cache.keys()
    for term in missing_terms:
        ancestor_cache[term] = go.get_ancestors(term)

    temp_files = []
    async_results = []

    model.eval()
    pool = Pool(processes=num_workers)

    try:
        with torch.no_grad():
            for chunk_id, chunk_start in enumerate(range(0, num_proteins, chunk_size)):
                chunk_end = min(chunk_start + chunk_size, num_proteins)
                chunk_protein_ids = protein_ids[chunk_start:chunk_end]

                chunk_predictions = []

                for batch_start in range(0, len(chunk_protein_ids), batch_size):
                    batch_end = min(batch_start + batch_size, len(chunk_protein_ids))
                    batch_protein_ids = chunk_protein_ids[batch_start:batch_end]

                    batch_features = []
                    for protein_id in batch_protein_ids:
                        vec = embeddings[protein_id]
                        if vec.shape != (vec_dim,):
                            raise ValueError(f"Embedding dimension differs for {protein_id}: {vec.shape}")
                        batch_features.append(vec)

                    batch_features_np = np.vstack(batch_features)
                    X_batch = torch.from_numpy(batch_features_np).to(device)
                    preds = model(X_batch).detach().cpu().numpy()

                    for i, protein_id in enumerate(batch_protein_ids):
                        if sparsity_threshold > 0:
                            mask = preds[i] >= sparsity_threshold
                            scores = {terms[j]: float(preds[i, j])
                                      for j in np.where(mask)[0]}
                        else:
                            scores = {terms[j]: float(preds[i, j])
                                      for j in range(len(terms))}
                        chunk_predictions.append((protein_id, scores))

                    del X_batch, preds, batch_features, batch_features_np

                temp_file = os.path.join(os.path.dirname(output_file),
                                         f".temp_{taxon_id}_{chunk_id:04d}.tsv")
                temp_files.append(temp_file)

                async_result = pool.apply_async(
                    _process_chunk_and_write,
                    ((chunk_id, chunk_predictions, ancestor_cache, sparsity_threshold, temp_file),)
                )
                async_results.append(async_result)
                if len(async_results) >= 2 * num_workers:
                    async_results.pop(0).get()

                del chunk_predictions

        for result in async_results:
            result.get()

    finally:
        pool.close()
        pool.join()

    partial = output_file + '.part'
    with open(partial, 'w') as out:
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                with open(temp_file, 'r') as inp:
                    shutil.copyfileobj(inp, out)
                os.remove(temp_file)

    os.replace(partial, output_file)
    del embeddings

    return num_proteins


@ck.command()
@ck.option('--subontology', help='Subontology.', default='cc', type=ck.Choice(['mf', 'bp', 'cc']))
@ck.option('--train-file', help='Path to train_proteins.tsv.',
           default='splits/heldout/train_proteins.tsv')
@ck.option('--val-file', help='Path to val_proteins.tsv.',
           default='splits/heldout/val_proteins.tsv')
@ck.option('--test-organisms-file', help='Path to test_organisms.txt.',
           default='splits/heldout/test_organisms.txt')
@ck.option('--esm2-embeddings-dir', help='Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl files.',
           default='${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings')
@ck.option('--output-dir', help='Output directory.', default='mlp_heldout_results')
@ck.option('--go-obo', help='GO OBO file.', default='data/go-basic.obo')
@ck.option('--terms-file', type=ck.Path(exists=True), default=None, help='Canonical terms.pkl for this GO aspect.')
@ck.option('--device', type=str, default='cuda:0' if torch.cuda.is_available() else 'cpu', help='Device')
@ck.option('--vec-dim', type=int, default=2560, help='Vector dimension')
@ck.option('--epochs', type=int, default=10, help='Training epochs')
@ck.option('--batch-size', type=int, default=256, help='Batch size')
@ck.option('--lr', type=float, default=1e-3, help='Learning rate')
@ck.option('--sparsity-threshold', type=float, default=0.0, help='Minimum score to include in output TSVs')
@ck.option('--num-workers', type=int, default=4, help='Number of workers for parallel processing')
@ck.option('--no-save-predictions', is_flag=True, default=False, help='Retired option; test evaluation requires exported predictions.')
@ck.option('--annotations-go-basic-dir',
           help='Retired option; supply source-rebuilt references to the shared evaluator.',
           default=None)
@ck.option('--fold-id', type=int, default=1, help='Fold ID for output filenames (predictions_fold_XX_taxon_YYYY.tsv). Compatible with run_adjustment_pipeline.py.')
@ck.option('--seed', type=int, default=42, help='Random seed for reproducible training (Python/NumPy/PyTorch).')
def main(subontology, train_file, val_file, test_organisms_file, esm2_embeddings_dir,
         output_dir, go_obo, terms_file, device, vec_dim, epochs, batch_size, lr,
         sparsity_threshold, num_workers, no_save_predictions, annotations_go_basic_dir, fold_id, seed):
    """
    Train MLP using a fixed train/val/test split and generate predictions for test organisms.
    """
    if no_save_predictions or annotations_go_basic_dir is not None:
        raise ck.UsageError(
            "Legacy annotated-only evaluation is retired. Export all test predictions "
            "and use analysis/aggregate_revision_audit.py with source-rebuilt references."
        )
    EMBEDDING_INPUTS.clear()
    set_seed(seed)
    print(f"[INFO] RNG seed set to {seed} for reproducible training")
    # Resolve absolute paths
    train_file = os.path.abspath(train_file)
    val_file = os.path.abspath(val_file)
    test_organisms_file = os.path.abspath(test_organisms_file)
    esm2_embeddings_dir = os.path.abspath(esm2_embeddings_dir)
    output_dir = os.path.abspath(output_dir)
    go_obo = os.path.abspath(go_obo)

    # Validate inputs
    for path, name in [(train_file, 'train_file'), (val_file, 'val_file'),
                       (test_organisms_file, 'test_organisms_file'), (go_obo, 'go_obo')]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{name} not found: {path}")

    os.makedirs(output_dir, exist_ok=True)

    # Load GO terms
    terms_file = terms_file or f'data/{subontology}/terms.pkl'
    if not os.path.exists(terms_file):
        raise FileNotFoundError(f"Terms file not found: {terms_file}")

    all_terms_array = pd.read_pickle(terms_file)['gos'].values.flatten()
    all_terms = list(all_terms_array)
    ck.echo(f"[INFO] Loaded {len(all_terms)} GO terms for {subontology}")

    training_manifest_path = os.path.join(os.path.dirname(train_file), 'training-manifest.json')
    if not os.path.isfile(training_manifest_path):
        raise FileNotFoundError('Run data_prep/prepare_baseline_inputs.py first; training-manifest.json is required')
    with open(training_manifest_path) as handle:
        training_manifest = json.load(handle)
    manifest_root = os.path.dirname(training_manifest_path)
    if training_manifest['go_sha256'] != digest(go_obo):
        raise ValueError("Training manifest uses another ontology")
    for path in [train_file, val_file, terms_file]:
        key = os.path.relpath(os.path.abspath(path), manifest_root)
        if training_manifest['files'].get(key) != digest(path):
            raise ValueError(f"Unverified corrected training input: {path}")

    # Load GO ontology
    go = Ontology(go_obo, with_rels=True)
    ck.echo(f"[INFO] Loaded GO ontology from {go_obo}")
    if len(all_terms) != len(set(all_terms)) or any(go.canonical_id(t) != t for t in all_terms):
        raise ValueError("Training vocabulary must contain unique, current primary GO IDs")

    # Load train and val DataFrames
    ck.echo(f"[INFO] Loading training proteins from {train_file}...")
    train_df = load_proteins_from_split_tsv(train_file)
    train_df = train_df[train_df['prop_annotations'].apply(len) > 0].reset_index(drop=True)
    ck.echo(f"[INFO] Loaded {len(train_df)} training proteins with annotations")

    ck.echo(f"[INFO] Loading validation proteins from {val_file}...")
    val_df = load_proteins_from_split_tsv(val_file)
    val_df = val_df[val_df['prop_annotations'].apply(len) > 0].reset_index(drop=True)
    ck.echo(f"[INFO] Loaded {len(val_df)} validation proteins with annotations")

    train_ids, val_ids = set(train_df.proteins), set(val_df.proteins)
    with open(test_organisms_file) as handle:
        test_taxa = set(handle.read().split())
    if (len(train_ids) != len(train_df) or len(val_ids) != len(val_df)
            or train_ids & val_ids or (set(train_df.orgs) | set(val_df.orgs)) & test_taxa):
        raise ValueError("Duplicate proteins or train/validation/test overlap")
    if train_df.empty or val_df.empty:
        raise ValueError("Empty train/validation split")

    # Train MLP
    ck.echo("[INFO] Training MLP...")
    start_time = time.time()
    model, train_info = train_mlp(train_df, all_terms, device, esm2_embeddings_dir, vec_dim,
                                   epochs=epochs, batch_size=batch_size, lr=lr)
    ck.echo(f"[INFO] Training completed in {time.time() - start_time:.2f} seconds")

    # Save the fitted state before potentially long validation and inference.
    torch.save({'state_dict': model.state_dict(), 'terms': all_terms,
                'vec_dim': vec_dim, 'seed': seed, 'epochs': epochs,
                'batch_size': batch_size, 'learning_rate': lr},
               os.path.join(output_dir, 'model.pt'))

    # Find optimal threshold on validation set
    ck.echo("[INFO] Finding optimal threshold on validation set...")
    val_metrics = find_optimal_threshold_on_validation(model, val_df, all_terms,
                                                       device, esm2_embeddings_dir, vec_dim)
    optimal_threshold = val_metrics['optimal_threshold']

    # Load test organism IDs
    test_taxon_ids: List[str] = []
    with open(test_organisms_file, 'r') as f:
        for line in f:
            taxon_id = line.strip()
            if taxon_id and not taxon_id.startswith('#'):
                test_taxon_ids.append(taxon_id)
    ck.echo(f"[INFO] Loaded {len(test_taxon_ids)} test organisms from {test_organisms_file}")

    # Generate predictions for each test organism
    predictions_dir = os.path.join(output_dir, 'predictions')
    os.makedirs(predictions_dir, exist_ok=True)

    ancestor_cache: Dict[str, Set[str]] = {}
    total_test_proteins = 0
    start_time = time.time()

    for taxon_idx, taxon_id in enumerate(test_taxon_ids, 1):
        ck.echo(f"[INFO] Processing taxon {taxon_idx}/{len(test_taxon_ids)}: {taxon_id}")

        output_file = os.path.join(predictions_dir, f"predictions_fold_{fold_id:02d}_taxon_{taxon_id}.tsv")
        taxon_start = time.time()

        num_proteins = generate_propagate_write_parallel(
            model=model,
            esm2_embeddings_dir=esm2_embeddings_dir,
            taxon_id=taxon_id,
            terms=all_terms,
            device=device,
            vec_dim=vec_dim,
            go=go,
            ancestor_cache=ancestor_cache,
            output_file=output_file,
            sparsity_threshold=sparsity_threshold,
            batch_size=batch_size,
            num_workers=num_workers,
            chunk_size=5000,
        )

        if num_proteins == 0:
            raise ValueError(f"No predictions generated for taxon {taxon_id}")

        ck.echo(f"[INFO] Processed {num_proteins} proteins for taxon {taxon_id} "
                f"in {time.time() - taxon_start:.2f} seconds")
        total_test_proteins += num_proteins

        if taxon_idx % 10 == 0:
            gc.collect()

    ck.echo(f"[INFO] Generated predictions for {total_test_proteins} proteins across "
            f"{len(test_taxon_ids)} test organisms in {time.time() - start_time:.2f} seconds")

    # Save results summary
    results = {
        'subontology': subontology,
        'n_train_proteins': len(train_df),
        'n_val_proteins': len(val_df),
        'n_test_proteins': total_test_proteins,
        'n_test_organisms': len(test_taxon_ids),
        'n_terms': len(all_terms),
        'validation_metric_definition': 'micro-averaged diagnostics on dense model scores; paper metrics are computed separately',
        'training_manifest_sha256': digest(training_manifest_path),
        'embedding_inputs': dict(sorted(EMBEDDING_INPUTS.items())),
        'software': {'python':platform.python_version(),'torch':torch.__version__,
                     'cuda':torch.version.cuda,'numpy':np.__version__,'pandas':pd.__version__,
                     'device':str(device),'gpu':torch.cuda.get_device_name(device) if str(device).startswith('cuda') else None},
        'code_sha256': {p:digest(os.path.join(os.path.dirname(__file__),'..','..',p)) for p in
                       ['baselines/mlp_esm2/train_mlp_heldout.py','deepgo/base.py','gaef_patches/ontology.py','data_prep/prepare_baseline_inputs.py']},
        'seed': seed, 'epochs': epochs, 'batch_size': batch_size, 'learning_rate': lr,
        'sparsity_threshold': sparsity_threshold,
        'checkpoint_sha256': digest(os.path.join(output_dir, 'model.pt')),
        **{f'val_{k}': v for k, v in val_metrics.items()},
        **{f'train_{k}': v for k, v in train_info.items()},
    }

    results_json = os.path.join(output_dir, 'results.json')
    with open(results_json, 'w') as f:
        json.dump(results, f, indent=2)

    results_csv = os.path.join(output_dir, 'results.csv')
    pd.DataFrame([results]).to_csv(results_csv, index=False)

    ck.echo(f"[INFO] Saved results to {results_json} and {results_csv}")
    ck.echo(f"\n=== SUMMARY ===")
    ck.echo(f"Subontology:        {subontology}")
    ck.echo(f"Train proteins:     {len(train_df)}")
    ck.echo(f"Val proteins:       {len(val_df)}")
    ck.echo(f"Test proteins:      {total_test_proteins}")
    ck.echo(f"Test organisms:     {len(test_taxon_ids)}")
    ck.echo(f"Optimal threshold:  {optimal_threshold:.4f}")
    ck.echo(f"\n--- Validation metrics ---")
    ck.echo(f"Fmax:               {val_metrics['fmax']:.4f}")
    ck.echo(f"Precision at Fmax:  {val_metrics['precision']:.4f}")
    ck.echo(f"Recall at Fmax:     {val_metrics['recall']:.4f}")
    ck.echo(f"AUC (micro):        {val_metrics['auc']:.4f}")


if __name__ == "__main__":
    main()
