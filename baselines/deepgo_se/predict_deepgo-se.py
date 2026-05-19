"""
Script to predict with pre-trained DeepGOZero ESM+ models on UniProt proteomes. (DeepGO-SE)
Combines DeepGOModel ensemble architecture from deepgo2/predict.py with
UniProt data pipeline from predict_mlp.py.

Given:
- data_root: contains deepgozero_esm_plus_*.th models, terms.pkl, go-plus.norm
- esm2_embeddings_dir: contains esm2_embeddings/<taxon_id>/esm2.pkl files
- output_dir: to save predictions for each test organism

Output format: protein_id<tab>GO:term|score<tab>...
"""

import click as ck
import os
import json
from typing import Dict, List, Tuple, Set
from collections import defaultdict
import time
import numpy as np

# NumPy version compatibility shim
# Files pickled with NumPy 2.x use numpy._core, older versions use numpy.core
import sys
try:
    import numpy._core.numeric  # NumPy 2.x
except ModuleNotFoundError:
    # Running on older NumPy - create module aliases for pickle compatibility
    import numpy.core as _np_core
    sys.modules['numpy._core'] = _np_core
    sys.modules['numpy._core.numeric'] = _np_core.numeric
    sys.modules['numpy._core.multiarray'] = _np_core.multiarray
    sys.modules['numpy._core._multiarray_umath'] = getattr(_np_core, '_multiarray_umath', _np_core.multiarray)

import pandas as pd
import torch
import torch.nn as nn
from tqdm import tqdm
import gc
import shutil
from multiprocessing import Pool

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))  # to have access to deepgo modules

from deepgo.models import DeepGOModel
from deepgo.data import load_normal_forms
from deepgo.utils import Ontology

# Route click echo through tqdm to avoid breaking progress bars
try:
    _original_ck_echo = ck.echo
    def _tqdm_echo(message: str = "", file=None, nl: bool = True, err: bool = False, color=None):
        tqdm.write(str(message))
    ck.echo = _tqdm_echo  # type: ignore
except Exception:
    pass


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
        ck.echo(f"[WARNING] ESM2 file not found for taxon {taxon_id}: {esm2_file}")
        return {}
    
    try:
        df = pd.read_pickle(esm2_file)
        if not {'protein', 'esm2'}.issubset(set(df.columns)):
            ck.echo(f"[WARNING] ESM2 file for taxon {taxon_id} missing required columns (protein, esm2)")
            return {}
        
        embeddings = {}
        for row in df.itertuples():
            protein_id = str(row.protein)
            vec = np.asarray(row.esm2, dtype=np.float32)
            embeddings[protein_id] = vec
        
        return embeddings
    except Exception as e:
        ck.echo(f"[WARNING] Error loading ESM2 file for taxon {taxon_id}: {e}")
        return {}


def generate_predictions_deepgozero(models: List[nn.Module], esm2_embeddings_dir: str, 
                                    taxon_id: str, terms: List[str], device: str, 
                                    vec_dim: int, sparsity_threshold: float = 0.0, 
                                    batch_size: int = 256) -> Dict[str, Dict[str, float]]:
    """
    Generate predictions for all proteins in a taxon using DeepGOZero ensemble models.
    Averages predictions across all models in the ensemble.
    
    Args:
        models: List of trained DeepGOModel instances
        esm2_embeddings_dir: Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl files
        taxon_id: Taxon ID
        terms: List of GO terms
        device: Device for model inference
        vec_dim: Vector dimension (2560 for ESM2)
        sparsity_threshold: Minimum score to include in predictions
        batch_size: Number of proteins to process in each batch
        
    Returns:
        Dictionary mapping protein_id -> {GO_term: score}
    """
    predictions = {}
    
    # Load embeddings for this taxon
    embeddings = load_esm2_embeddings_taxon(esm2_embeddings_dir, taxon_id)
    if not embeddings:
        return predictions
    
    # Get protein IDs
    protein_ids = list(embeddings.keys())
    num_proteins = len(protein_ids)
    
    # Initialize sum of predictions for ensemble averaging
    sum_preds = np.zeros((num_proteins, len(terms)), dtype=np.float32)
    
    # Process each model in the ensemble
    for model_idx, model in enumerate(models):
        model.eval()
        
        with torch.no_grad():
            # Process proteins in batches
            for batch_start in tqdm(range(0, num_proteins, batch_size), 
                                   desc=f"Model {model_idx+1}/{len(models)} - Taxon {taxon_id}", 
                                   total=(num_proteins + batch_size - 1) // batch_size):
                batch_end = min(batch_start + batch_size, num_proteins)
                batch_protein_ids = protein_ids[batch_start:batch_end]
                
                # Prepare features for this batch
                batch_features = []
                for protein_id in batch_protein_ids:
                    vec = embeddings[protein_id]
                    if vec.shape[0] == vec_dim:
                        batch_features.append(vec)
                    elif vec.shape[0] > vec_dim:
                        batch_features.append(vec[:vec_dim])
                    else:
                        padded = np.zeros(vec_dim, dtype=np.float32)
                        padded[:vec.shape[0]] = vec
                        batch_features.append(padded)
                
                # Stack features and run inference
                batch_features_np = np.vstack(batch_features)
                X_batch = torch.from_numpy(batch_features_np).to(device)
                
                # Run inference on batch
                logits = model(X_batch)
                preds = logits.detach().cpu().numpy()
                
                # Accumulate predictions
                sum_preds[batch_start:batch_end] += preds
                
                # Clean up batch tensors
                del X_batch, logits, preds
    
    # Average predictions across ensemble
    avg_preds = sum_preds / len(models)
    
    # Build predictions dictionary
    if sparsity_threshold > 0:
        # Apply threshold filter using vectorized operation
        mask = avg_preds >= sparsity_threshold
        
        for i, protein_id in enumerate(protein_ids):
            # Get indices where predictions exceed threshold
            term_indices = np.where(mask[i])[0]
            if len(term_indices) > 0:
                # Create dictionary for this protein
                predictions[protein_id] = {
                    terms[j]: float(avg_preds[i, j]) 
                    for j in term_indices
                }
    else:
        # No threshold - include all predictions
        for i, protein_id in enumerate(protein_ids):
            predictions[protein_id] = {
                terms[j]: float(avg_preds[i, j]) 
                for j in range(len(terms))
            }
    
    return predictions


def _propagate_single_protein(protein_id: str, go_scores: Dict[str, float],
                               ancestor_cache: Dict[str, Set[str]]) -> Dict[str, float]:
    """
    Propagate scores for a single protein through ontology hierarchy.

    Args:
        protein_id: Protein ID (unused but kept for consistency)
        go_scores: Dictionary mapping GO term to score
        ancestor_cache: Pre-computed ancestor cache

    Returns:
        Dictionary with propagated scores
    """
    propagated_scores: Dict[str, float] = {}

    for go_id, score in go_scores.items():
        # Update score for this term (take max)
        if go_id not in propagated_scores or propagated_scores[go_id] < score:
            propagated_scores[go_id] = score

        # Propagate to ancestors
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
            # Propagate scores through ontology
            propagated = _propagate_single_protein(protein_id, go_scores, ancestor_cache)

            # Format protein ID
            if '|' in protein_id:
                protein_id = protein_id.split('|')[2].split()[0]

            # Format and write line
            if not propagated:
                f.write(f"{protein_id}\n")
            else:
                filtered_terms = [f"{go_term}|{score:.6f}"
                                 for go_term, score in propagated.items()
                                 if score >= sparsity_threshold]
                if filtered_terms:
                    f.write(f"{protein_id}\t" + "\t".join(filtered_terms) + "\n")
                else:
                    f.write(f"{protein_id}\n")

    return temp_file


def generate_propagate_write_parallel(models: List[nn.Module], esm2_embeddings_dir: str, 
                                       taxon_id: str, terms: List[str], device: str, 
                                       vec_dim: int, go: Ontology, 
                                       ancestor_cache: Dict[str, Set[str]], output_file: str, 
                                       sparsity_threshold: float = 0.0, batch_size: int = 256, 
                                       num_workers: int = 4, chunk_size: int = 5000) -> int:
    """
    Generate predictions with DeepGOZero ensemble, propagate scores, and write to TSV 
    using parallel workers.

    Main process handles DeepGOZero inference (GPU-bound).
    Worker processes handle propagation and writing to temp files (CPU-bound).
    Temp files are merged at the end.

    Args:
        models: List of trained DeepGOModel instances
        esm2_embeddings_dir: Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl files
        taxon_id: Taxon ID
        terms: List of GO terms
        device: Device for model inference
        vec_dim: Vector dimension
        go: Ontology object for ancestor lookups
        ancestor_cache: Pre-computed ancestor cache (will be updated)
        output_file: Path to output TSV file
        sparsity_threshold: Minimum score to include in output
        batch_size: Batch size for model inference
        num_workers: Number of parallel workers
        chunk_size: Number of proteins per chunk for parallel processing

    Returns:
        Number of proteins processed
    """
    # Load embeddings for this taxon
    embeddings = load_esm2_embeddings_taxon(esm2_embeddings_dir, taxon_id)
    if not embeddings:
        return 0

    protein_ids = list(embeddings.keys())
    num_proteins = len(protein_ids)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Pre-fetch ancestors for all terms (needed by workers)
    missing_terms = set(terms) - ancestor_cache.keys()
    for term in missing_terms:
        ancestor_cache[term] = go.get_ancestors(term)

    temp_files = []
    async_results = []

    pool = Pool(processes=num_workers)

    try:
        # Initialize sum of predictions for ensemble averaging
        sum_preds = np.zeros((num_proteins, len(terms)), dtype=np.float32)
        
        # Process each model in the ensemble
        for model_idx, model in enumerate(models):
            ck.echo(f"[INFO] Processing model {model_idx+1}/{len(models)} for taxon {taxon_id}")
            model.eval()
            
            with torch.no_grad():
                for batch_start in range(0, num_proteins, batch_size):
                    batch_end = min(batch_start + batch_size, num_proteins)
                    batch_protein_ids = protein_ids[batch_start:batch_end]
                    
                    # Prepare features for this batch
                    batch_features = []
                    for protein_id in batch_protein_ids:
                        vec = embeddings[protein_id]
                        if vec.shape[0] == vec_dim:
                            batch_features.append(vec)
                        elif vec.shape[0] > vec_dim:
                            batch_features.append(vec[:vec_dim])
                        else:
                            padded = np.zeros(vec_dim, dtype=np.float32)
                            padded[:vec.shape[0]] = vec
                            batch_features.append(padded)
                    
                    # Run inference
                    batch_features_np = np.vstack(batch_features)
                    X_batch = torch.from_numpy(batch_features_np).to(device)
                    preds = model(X_batch).detach().cpu().numpy()
                    
                    # Accumulate predictions
                    sum_preds[batch_start:batch_end] += preds
                    
                    # Free batch memory
                    del X_batch, preds, batch_features, batch_features_np
        
        # Average predictions across ensemble
        avg_preds = sum_preds / len(models)
        
        # Process proteins in chunks with parallel workers
        for chunk_id, chunk_start in enumerate(range(0, num_proteins, chunk_size)):
            chunk_end = min(chunk_start + chunk_size, num_proteins)
            chunk_protein_ids = protein_ids[chunk_start:chunk_end]
            
            # Build predictions for this chunk (apply sparsity threshold early)
            chunk_predictions = []
            for i, protein_id in enumerate(chunk_protein_ids):
                idx = chunk_start + i
                if sparsity_threshold > 0:
                    # Only keep scores above threshold
                    mask = avg_preds[idx] >= sparsity_threshold
                    scores = {terms[j]: float(avg_preds[idx, j])
                              for j in np.where(mask)[0]}
                else:
                    scores = {terms[j]: float(avg_preds[idx, j])
                              for j in range(len(terms))}
                chunk_predictions.append((protein_id, scores))
            
            # Submit chunk to worker (non-blocking)
            temp_file = os.path.join(os.path.dirname(output_file),
                                    f".temp_{taxon_id}_{chunk_id:04d}.tsv")
            temp_files.append(temp_file)
            
            async_result = pool.apply_async(
                _process_chunk_and_write,
                ((chunk_id, chunk_predictions, ancestor_cache, sparsity_threshold, temp_file),)
            )
            async_results.append(async_result)
            
            # Free chunk memory immediately after dispatching
            del chunk_predictions

        # Wait for all workers to complete
        for result in async_results:
            result.get()  # This will raise any exceptions from workers

    finally:
        pool.close()
        pool.join()

    # Merge temp files in order
    with open(output_file, 'w') as out:
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                with open(temp_file, 'r') as inp:
                    shutil.copyfileobj(inp, out)
                os.remove(temp_file)

    # Clean up embeddings
    del embeddings

    return num_proteins


@ck.command()
@ck.option('--subontology', help='Subontology.', default='cc', type=ck.Choice(['mf', 'bp', 'cc']))
@ck.option('--data-root', help='Data root directory with GO OBO file.',
           default='data')
@ck.option('--models-dir', help='Models directory with DeepGOZero models.',
           default='deepgozero_models')
@ck.option('--esm2-embeddings-dir', help='Base directory containing esm2_embeddings/<taxon_id>/esm2.pkl files.',
           default='${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings')
@ck.option('--test-organisms-file', help='Text file with one test-organism taxon ID per line.',
           default='splits/timeset/test_organisms.txt')
@ck.option('--output-dir', help='Output directory.',
           default='deepgo-se_results')
@ck.option('--fold-id', type=int, default=1,
           help='Fold ID written into predictions_fold_XX_taxon_YYYY.tsv filenames. Kept for compatibility with pipeline/run_adjustment_pipeline.py.')
@ck.option('--device', type=str, default='cuda:0' if torch.cuda.is_available() else 'cpu', 
           help='Device')
@ck.option('--vec-dim', type=int, default=2560, help='Vector dimension (ESM2 default: 2560)')
@ck.option('--batch-size', type=int, default=256, help='Batch size for inference')
@ck.option('--sparsity-threshold', type=float, default=0.0, 
           help='Minimum score to include in output TSVs')
@ck.option('--go-obo', help='GO OBO file.', default='data/go-basic.obo')
@ck.option('--num-workers', type=int, default=4, 
           help='Number of workers for parallel processing')
def main(subontology, data_root, models_dir, esm2_embeddings_dir, test_organisms_file, output_dir, fold_id,
         device, vec_dim, batch_size, sparsity_threshold, go_obo, num_workers):
    """
    Main function for DeepGOZero ESM+ prediction on UniProt proteomes.
    """
    # Validate and normalize input paths
    try:
        data_root = os.path.abspath(data_root)
        models_dir = os.path.abspath(models_dir)
        esm2_embeddings_dir = os.path.abspath(esm2_embeddings_dir)
        test_organisms_file = os.path.abspath(test_organisms_file)
        output_dir = os.path.abspath(output_dir)

        if not os.path.exists(data_root):
            raise FileNotFoundError(f"Data root directory not found: {data_root}")
        if not os.path.exists(models_dir):
            raise FileNotFoundError(f"Models directory not found: {models_dir}")
        if not os.path.exists(test_organisms_file):
            raise FileNotFoundError(f"Test organisms file not found: {test_organisms_file}")
    except Exception as e:
        ck.echo(f"[ERROR] Error validating input paths: {e}")
        raise
    
    # Load terms
    terms_file = f'{data_root}/{subontology}/terms.pkl'
    if not os.path.exists(terms_file):
        raise FileNotFoundError(f"Terms file not found: {terms_file}")
    
    try:
        all_terms_array = pd.read_pickle(terms_file)['gos'].values.flatten()
        all_terms = list(all_terms_array)
        ck.echo(f"[INFO] Loaded {len(all_terms)} GO terms for {subontology}")
    except Exception as e:
        raise ValueError(f"[ERROR] Error loading terms from {terms_file}: {e}")
    
    # Create terms dictionary
    terms_dict = {v: i for i, v in enumerate(all_terms)}
    
    # Load GO normal forms
    go_norm = f'{data_root}/go-plus.norm'
    if not os.path.exists(go_norm):
        raise FileNotFoundError(f"GO normal forms file not found: {go_norm}")
    
    ck.echo(f"[INFO] Loading GO normal forms from {go_norm}")
    nf1, nf2, nf3, nf4, relations, zero_classes = load_normal_forms(go_norm, terms_dict)
    n_rels = len(relations)
    n_zeros = len(zero_classes)
    ck.echo(f"[INFO] Loaded normal forms: {len(nf1)} NF1, {len(nf2)} NF2, {len(nf3)} NF3, {len(nf4)} NF4")
    ck.echo(f"[INFO] Relations: {n_rels}, Zero-shot classes: {n_zeros}")
    
    # Define ensemble model indices per ontology
    ent_models = {
        'mf': [0, 1, 2, 5, 6, 8],
        'bp': [2, 5, 6, 7, 8, 9],
        'cc': [1, 3, 4, 5, 6, 7]
    }
    
    # Load ensemble models
    models = []
    for mn in ent_models[subontology]:
        model_file = f'{models_dir}/deepgozero_esm_plus_{mn}.th'
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        
        model = DeepGOModel(vec_dim, len(all_terms), n_zeros, n_rels, device).to(device)
        model.load_state_dict(torch.load(model_file, map_location=device))
        model.eval()
        models.append(model)
        ck.echo(f"[INFO] Loaded model {mn}: {model_file}")
    
    ck.echo(f"[INFO] Loaded {len(models)} ensemble models for {subontology}")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        raise OSError(f"[ERROR] Failed to create output directory {output_dir}: {e}")

    with open(test_organisms_file, 'r') as f:
        test_organisms = [line.strip() for line in f if line.strip()]
    if not test_organisms:
        raise ValueError(f"No taxon IDs read from {test_organisms_file}")
    ck.echo(f"[INFO] Loaded {len(test_organisms)} test organisms from {test_organisms_file}")

    # Load GO ontology for score propagation
    if not os.path.exists(go_obo):
        raise FileNotFoundError(f"GO OBO file not found: {go_obo}")
    go = Ontology(go_obo, with_rels=True)
    ck.echo(f"[INFO] Loaded GO ontology from {go_obo}")

    target_dir = os.path.join(output_dir, 'predictions')
    os.makedirs(target_dir, exist_ok=True)

    ck.echo(f"[INFO] Generating predictions for {len(test_organisms)} test organisms...")
    start_time = time.time()
    ancestor_cache: Dict = {}
    total_test_predictions = 0
    written_taxons: List[str] = []

    for taxon_idx, taxon_id in enumerate(test_organisms, 1):
        ck.echo(f"[INFO] Processing taxon {taxon_idx}/{len(test_organisms)}: {taxon_id}")
        output_file = os.path.join(target_dir, f"predictions_fold_{fold_id:02d}_taxon_{taxon_id}.tsv")
        taxon_start_time = time.time()

        num_proteins = generate_propagate_write_parallel(
            models=models,
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
            ck.echo(f"[WARNING] No predictions generated for taxon {taxon_id}")
            continue

        ck.echo(f"[INFO] Processed {num_proteins} proteins for taxon {taxon_id} in {time.time() - taxon_start_time:.2f} seconds")
        total_test_predictions += num_proteins
        written_taxons.append(taxon_id)

        if taxon_idx % 10 == 0:
            gc.collect()

    ck.echo(f"[INFO] Generated predictions for {total_test_predictions} proteins from {len(written_taxons)} test organisms in {time.time() - start_time:.2f} seconds")

    summary = {
        'subontology': subontology,
        'fold_id': fold_id,
        'test_organisms_file': test_organisms_file,
        'n_test_organisms_requested': len(test_organisms),
        'n_test_organisms_with_preds': len(written_taxons),
        'n_test_predictions': total_test_predictions,
    }
    results_json = os.path.join(output_dir, f'cv_results_fold_{fold_id:02d}.json')
    with open(results_json, 'w') as f:
        json.dump(summary, f, indent=2)

    ck.echo("[INFO] Done. Saved test organism predictions.")
    ck.echo(f"  - Test Organisms: {len(written_taxons)}")
    ck.echo(f"  - Test Proteins: {total_test_predictions}")


if __name__ == "__main__":
    main()
