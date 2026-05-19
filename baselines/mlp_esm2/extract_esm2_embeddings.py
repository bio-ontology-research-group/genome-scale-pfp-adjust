#!/usr/bin/env python3
"""
Script to extract ESM2 embeddings for all proteins in the UniProt proteomes.

Given:
- uniprot_proteomes_ids.tsv file
- uniprot_proteomes_dir directory:
    - contains the protein sequences in the format: INPUT_DIR/<domain>/<proteome_id>/<proteome_id>_<taxon_id>.fasta.gz

Output:
- esm2.pkl file for each taxon saved under OUTPUT_DIR/esm2_embeddings/<taxon_id>/ directory
"""
import os
import sys
import time
import signal
from typing import List, Tuple, Optional
import gzip
import psutil
import click as ck
import pandas as pd
import torch
import numpy as np
from esm import FastaBatchedDataset, pretrained

# based on deepgo/extract_esm.py
class GzippedFastaBatchedDataset(FastaBatchedDataset):
    @classmethod
    def from_file(cls, fasta_file):
        sequence_labels, sequence_strs = [], []
        cur_seq_label = None
        buf = []

        def _flush_current_seq():
            nonlocal cur_seq_label, buf
            if cur_seq_label is None:
                return
            sequence_labels.append(cur_seq_label)
            sequence_strs.append("".join(buf))
            cur_seq_label = None
            buf = []

        with gzip.open(fasta_file, "rt") as infile:
            for line_idx, line in enumerate(infile):
                if line.startswith(">"):  # label line
                    _flush_current_seq()
                    line = line[1:].strip()
                    if len(line) > 0:
                        cur_seq_label = line
                    else:
                        cur_seq_label = f"seqnum{line_idx:09d}"
                else:  # sequence line
                    buf.append(line.strip())

        _flush_current_seq()

        assert len(set(sequence_labels)) == len(
            sequence_labels
        ), "Found duplicate sequence labels"

        return cls(sequence_labels, sequence_strs)


def load_proteomes_metadata(tsv_file: str) -> List[Tuple[str, str, str]]:
    """
    Load and parse the TSV file, return list of (taxon_id, proteome_id, domain) tuples.
    
    Args:
        tsv_file: Path to uniprot_proteomes_ids.tsv
        
    Returns:
        List of (taxon_id, proteome_id, domain) tuples
    """
    df = pd.read_csv(tsv_file, sep='\t')
    metadata = []
    for _, row in df.iterrows():
        taxon_id = str(row['TaxonID'])
        proteome_id = str(row['ProteomeID'])
        domain = str(row['Domain']).lower()
        metadata.append((taxon_id, proteome_id, domain))
    return metadata


def find_proteome_fasta_file(input_dir: str, domain: str, proteome_id: str, taxon_id: str) -> Optional[str]:
    """
    Locate the FASTA file for a given proteome.
    Handles domain capitalization: 'bacteria' -> 'Bacteria', etc.
    
    Args:
        input_dir: Base directory for proteome FASTA files
        domain: Domain name (lowercase from TSV)
        proteome_id: UniProt proteome ID
        taxon_id: Taxon ID
        
    Returns:
        Path to FASTA file or None if not found
    """
    # Map lowercase domain to capitalized directory name
    domain_map = {
        'bacteria': 'Bacteria',
        'eukaryota': 'Eukaryota',
        'archaea': 'Archaea',
        'viruses': 'Viruses'
    }
    
    capitalized_domain = domain_map.get(domain.lower(), domain.capitalize())
    
    # Try the expected path
    fasta_file = os.path.join(
        input_dir, 
        capitalized_domain, 
        proteome_id, 
        f'{proteome_id}_{taxon_id}.fasta.gz'
    )
    
    if os.path.exists(fasta_file):
        return fasta_file
    
    # Fallback: try all domains if exact match fails
    for fallback_domain in ['Bacteria', 'Eukaryota', 'Archaea', 'Viruses']:
        fallback_file = os.path.join(
            input_dir,
            fallback_domain,
            proteome_id,
            f'{proteome_id}_{taxon_id}.fasta.gz'
        )
        if os.path.exists(fallback_file):
            return fallback_file
    
    return None


def load_esm_model(model_location='esm2_t36_3B_UR50D', device=None):
    """Load ESM2 model with proper error handling and logging."""
    print(f"Loading ESM2 model: {model_location}")
    
    start_time = time.time()
    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    load_time = time.time() - start_time
    
    print(f"Model loaded successfully in {load_time:.1f}s")
    
    model.eval()
    if device and device.startswith('cuda'):
        print(f"Moving model to device: {device}")
        model = model.to(device)
    
    return model, alphabet


def extract_esm2(fasta_file, model, alphabet,
                 truncation_seq_length=1022, toks_per_batch=4096,
                 device=None):
    """
    Extract ESM2 embeddings from a FASTA file.
    Uses GzippedFastaBatchedDataset for .fasta.gz files.
    
    Returns:
        Tuple of (proteins list, vectors numpy array)
    """
    # Use GzippedFastaBatchedDataset for .fasta.gz files
    if fasta_file.endswith('.gz'):
        dataset = GzippedFastaBatchedDataset.from_file(fasta_file)
    else:
        dataset = FastaBatchedDataset.from_file(fasta_file)
    
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )
    
    print(f"Processing {fasta_file} with {len(dataset)} sequences in {len(batches)} batches")
    
    proteins = []
    vectors = []
    repr_layers = [36,]
    
    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(f"Processing batch {batch_idx + 1}/{len(batches)} ({toks.size(0)} sequences)")
            
            if device:
                toks = toks.to(device, non_blocking=True)
            
            out = model(toks, repr_layers=repr_layers, return_contacts=False)
            
            # Get representations from layer 36
            final_repr = out["representations"][36]  # [batch_size, seq_len, hidden_dim]
            batch_size, seq_len, hidden_dim = final_repr.shape
            
            # Process each sequence
            for i in range(batch_size):
                # Get actual sequence length
                truncate_len = min(truncation_seq_length, len(strs[i]))
                
                # Extract mean representation (skip [CLS] token)
                mean_repr = final_repr[i, 1:truncate_len + 1].mean(0)  # [hidden_dim]
                
                # Move to CPU and convert to numpy
                if device and device != 'cpu':
                    mean_repr = mean_repr.cpu()
                
                proteins.append(labels[i])
                vectors.append(mean_repr.numpy())
            
            
            # Clear GPU cache periodically
            if device and device != 'cpu' and (batch_idx + 1) % 10 == 0:
                torch.cuda.empty_cache()
    
    # Stack all vectors into a single numpy array
    vectors_array = np.vstack(vectors).astype(np.float32)
    
    print(f"Completed processing {len(proteins)} proteins from {fasta_file}")
    return proteins, vectors_array


def process_taxon(taxid: str, proteome_id: str, domain: str, input_dir: str, output_dir: str, 
                 output_name: str, device: str, overwrite: bool, batch_size: int, 
                 model=None, alphabet=None) -> Tuple[str, bool, int, Optional[str]]:
    """
    Process a single taxon and save ESM2 embeddings.
    
    Returns:
        (taxon_id, is_skipped, protein_count, output_path_or_error)
    """
    try:
        # Find the FASTA file
        fasta_file = find_proteome_fasta_file(input_dir, domain, proteome_id, taxid)
        if not fasta_file:
            print(f"Warning: FASTA file not found for taxon {taxid} (proteome {proteome_id})")
            return taxid, True, 0, None  # skipped - no file
        
        # Create output directory for this taxon
        taxon_output_dir = os.path.join(output_dir, 'esm2_embeddings', taxid)
        os.makedirs(taxon_output_dir, exist_ok=True)
        out_path = os.path.join(taxon_output_dir, output_name)
        
        # Check if output exists and skip if not overwriting
        if (not overwrite) and os.path.exists(out_path):
            print(f"Skipping {taxid}: output already exists")
            return taxid, True, 0, None  # skipped
        
        print(f"Processing taxon {taxid} (proteome {proteome_id}, domain {domain})")
        print(f"FASTA file: {fasta_file}")
        
        # Extract ESM2 embeddings
        proteins, vectors = extract_esm2(
            fasta_file=fasta_file,
            model=model,
            alphabet=alphabet,
            device=device,
            toks_per_batch=batch_size
        )
        
        # Create DataFrame with 'protein' and 'esm2' columns
        pd.DataFrame({
            'protein': proteins,
            'esm2': [vec for vec in vectors]
        }).to_pickle(out_path)

        print(f"Saved {len(proteins)} proteins to {out_path}")
        # print memory usage
        print(f"Memory usage: {psutil.virtual_memory().percent}%")
        return taxid, False, len(proteins), out_path  # success
    except Exception as e:
        error_msg = f"Error processing taxon {taxid}: {str(e)}"
        print(error_msg)
        return taxid, False, 0, error_msg  # error


def setup_signal_handlers():
    """Setup signal handlers for graceful shutdown."""
    def signal_handler(signum, frame):
        print(f"\nReceived signal {signum}. Gracefully shutting down...")
        sys.exit(0)
    
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)


def save_checkpoint(output_dir: str, processed_taxids: List[str], start_time: float, 
                   processed_count: int, skipped_count: int, error_count: int, 
                   device: str = "unknown", run_id: str = None):
    """Save checkpoint information."""
    # Extract GPU ID from device string for unique checkpoint files
    gpu_id = device.split(':')[-1] if ':' in device else 'unknown'
    
    # Create unique checkpoint filename
    if run_id:
        checkpoint_file = os.path.join(output_dir, f"esm2_extraction_checkpoint_gpu_{gpu_id}_{run_id}.txt")
    else:
        checkpoint_file = os.path.join(output_dir, f"esm2_extraction_checkpoint_gpu_{gpu_id}.txt")
    
    checkpoint_data = {
        'timestamp': time.time(),
        'processed_taxids': processed_taxids,
        'start_time': start_time,
        'processed_count': processed_count,
        'skipped_count': skipped_count,
        'error_count': error_count,
        'elapsed_time': time.time() - start_time
    }
    
    with open(checkpoint_file, 'w') as f:
        for key, value in checkpoint_data.items():
            if key == 'processed_taxids':
                f.write(f"{key}: {','.join(map(str, value))}\n")
            else:
                f.write(f"{key}: {value}\n")


def load_checkpoint(output_dir: str, device: str = "unknown", run_id: str = None) -> Tuple[List[str], float, int, int, int]:
    """Load checkpoint information if available."""
    # Extract GPU ID from device string for unique checkpoint files
    gpu_id = device.split(':')[-1] if ':' in device else 'unknown'
    
    # Create unique checkpoint filename
    if run_id:
        checkpoint_file = os.path.join(output_dir, f"esm2_extraction_checkpoint_gpu_{gpu_id}_{run_id}.txt")
    else:
        checkpoint_file = os.path.join(output_dir, f"esm2_extraction_checkpoint_gpu_{gpu_id}.txt")
    
    if not os.path.exists(checkpoint_file):
        return [], time.time(), 0, 0, 0
    
    try:
        processed_taxids = []
        start_time = time.time()
        processed_count = 0
        skipped_count = 0
        error_count = 0
        
        with open(checkpoint_file, 'r') as f:
            for line in f:
                if ':' not in line:
                    continue
                key, value = line.strip().split(': ', 1)
                if key == 'processed_taxids':
                    processed_taxids = [tid.strip() for tid in value.split(',') if tid.strip()]
                elif key == 'start_time':
                    start_time = float(value)
                elif key == 'processed_count':
                    processed_count = int(value)
                elif key == 'skipped_count':
                    skipped_count = int(value)
                elif key == 'error_count':
                    error_count = int(value)
        
        print(f"Loaded checkpoint: {len(processed_taxids)} taxa already processed")
        return processed_taxids, start_time, processed_count, skipped_count, error_count
        
    except Exception as e:
        print(f"Warning: Could not load checkpoint: {e}")
        return [], time.time(), 0, 0, 0


def parse_taxon_ids(taxon_ids_arg: Optional[str]) -> Optional[List[str]]:
    """Parse taxon IDs from comma-separated string or file path."""
    if not taxon_ids_arg:
        return None
    
    # Check if it's a file path
    if os.path.exists(taxon_ids_arg):
        with open(taxon_ids_arg, 'r') as f:
            taxon_ids = [line.strip() for line in f if line.strip()]
        return taxon_ids
    
    # Otherwise, treat as comma-separated string
    taxon_ids = [tid.strip() for tid in taxon_ids_arg.split(',') if tid.strip()]
    return taxon_ids if taxon_ids else None


@ck.command()
@ck.option('--proteomes-ids-file', '-p', default='data/uniprot_proteomes_ids.tsv',
           help='Path to uniprot_proteomes_ids.tsv file.')
@ck.option('--input-dir', '-i', default='${DATA_DIR}/uniprot_reference_proteomes',
           help='Base directory containing <domain>/<proteome_id>/<proteome_id>_<taxon_id>.fasta.gz files.')
@ck.option('--output-dir', '-o', default='${DATA_DIR}/swissprot_proteomes_folds/esm2_embeddings',
           help='Base output directory where esm2_embeddings/<taxon_id>/esm2.pkl will be saved.')
@ck.option('--start-idx', type=int, default=None,
           help='Starting index in TSV (0-based, for parallel jobs).')
@ck.option('--end-idx', type=int, default=None,
           help='Ending index in TSV (exclusive, for parallel jobs).')
@ck.option('--taxon-ids', type=str, default=None,
           help='Comma-separated taxon IDs or file path (alternative to idx range).')
@ck.option('--device', '-d', default='cuda:0',
           help='Device for ESM2 model, e.g., cuda:0 or cpu')
@ck.option('--batch-size', '-bs', default=4096,
           help='Batch size (tokens per batch) for ESM2 processing.')
@ck.option('--overwrite', is_flag=True, default=False,
           help='Recompute even if output exists.')
@ck.option('--run-id', default=None,
           help='Unique identifier for this run to avoid file conflicts.')
@ck.option('--checkpoint-interval', default=1,
           help='Save checkpoint every N taxa processed.')
@ck.option('--output-name', default='esm2.pkl',
           help='Output filename (default: esm2.pkl).')
def main(proteomes_ids_file: str, input_dir: str, output_dir: str,
         start_idx: Optional[int], end_idx: Optional[int], taxon_ids: Optional[str],
         device: str, batch_size: int, overwrite: bool, run_id: str,
         checkpoint_interval: int, output_name: str):
    
    # Setup signal handlers for graceful shutdown
    setup_signal_handlers()
    
    # Check device availability
    if device.startswith('cuda') and not torch.cuda.is_available():
        print(f"Error: CUDA device {device} not available")
        return
    
    print(f"Processing on device: {device}")
    
    # Load ESM2 model and alphabet once
    model, alphabet = load_esm_model(device=device)
    
    # Load proteomes metadata
    all_metadata = load_proteomes_metadata(proteomes_ids_file)
    print(f"Loaded {len(all_metadata)} taxa from {proteomes_ids_file}")
    
    # Filter metadata based on start_idx/end_idx or taxon_ids
    if taxon_ids:
        # Process specific taxon IDs
        requested_taxids = parse_taxon_ids(taxon_ids)
        if not requested_taxids:
            print("Error: No valid taxon IDs provided")
            return
        
        # Create a set for fast lookup
        requested_set = set(requested_taxids)
        metadata = [(tid, pid, dom) for tid, pid, dom in all_metadata if tid in requested_set]
        print(f"Filtered to {len(metadata)} taxa from requested list")
    elif start_idx is not None and end_idx is not None:
        # Process range
        metadata = all_metadata[start_idx:end_idx]
        print(f"Processing range [{start_idx}:{end_idx}]: {len(metadata)} taxa")
    else:
        # Process all
        metadata = all_metadata
        print(f"Processing all {len(metadata)} taxa")
    
    if not metadata:
        print("No taxa to process")
        return
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load checkpoint if available
    processed_taxids, start_time, processed_count, skipped_count, error_count = load_checkpoint(
        output_dir, device, run_id
    )
    
    # Filter out already processed taxids
    if not overwrite and processed_taxids:
        original_count = len(metadata)
        processed_set = set(processed_taxids)
        metadata = [(tid, pid, dom) for tid, pid, dom in metadata if tid not in processed_set]
        print(f"Resuming from checkpoint: {len(metadata)} taxa remaining out of {original_count}")
    
    processed = processed_count
    skipped = skipped_count
    errors = error_count
    total_proteins = 0
    
    if not start_time:
        start_time = time.time()
    
    try:
        for i, (taxid, proteome_id, domain) in enumerate(metadata):
            print(f"\nProcessing taxon {i+1}/{len(metadata)}: {taxid}")
            
            result = process_taxon(
                taxid, proteome_id, domain, input_dir, output_dir,
                output_name, device, overwrite, batch_size,
                model=model, alphabet=alphabet
            )
            
            taxid_result, is_skipped, protein_count, output_path = result
            
            if is_skipped:
                if protein_count == 0:
                    skipped += 1
                else:
                    processed += 1
                    total_proteins += protein_count
            elif output_path is None or not isinstance(output_path, str) or not output_path.endswith(output_name):
                errors += 1
                print(f"Error processing {taxid_result}: {output_path}")
            else:
                processed += 1
                total_proteins += protein_count
            
            processed_taxids.append(taxid)
            
            # Save checkpoint periodically
            if (i + 1) % checkpoint_interval == 0:
                save_checkpoint(output_dir, processed_taxids, start_time, processed, skipped, errors, device, run_id)
                print(f"Checkpoint saved: {processed} processed, {skipped} skipped, {errors} errors")
            
            # Clear GPU cache after each taxon
            if device != 'cpu':
                torch.cuda.empty_cache()
                
    except KeyboardInterrupt:
        print("\nReceived interrupt signal. Saving checkpoint and exiting...")
        save_checkpoint(output_dir, processed_taxids, start_time, processed, skipped, errors, device, run_id)
        sys.exit(1)
    
    elapsed_time = time.time() - start_time
    
    print("\n" + "="*60)
    print("PROCESSING COMPLETE")
    print("="*60)
    print(f"Processed:  {processed} taxa")
    print(f"Skipped:    {skipped} taxa")
    print(f"Errors:     {errors} taxa")
    print(f"Total proteins: {total_proteins:,}")
    print(f"Elapsed time: {elapsed_time:.1f}s")
    if processed > 0:
        print(f"Avg time per taxon: {elapsed_time/processed:.1f}s")
        print(f"Avg time per protein: {elapsed_time/total_proteins*1000:.2f}ms")
    
    # Save summary
    device_clean = device.replace(':', '_')
    if run_id:
        summary_file = os.path.join(output_dir, f"esm2_extraction_summary_{device_clean}_{run_id}.tsv")
    else:
        summary_file = os.path.join(output_dir, f"esm2_extraction_summary_{device_clean}.tsv")
    
    summary_df = pd.DataFrame({
        'device': [device],
        'processed': [processed],
        'skipped': [skipped],
        'errors': [errors],
        'total_proteins': [total_proteins],
        'elapsed_time': [elapsed_time],
        'batch_size': [batch_size]
    })
    summary_df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary saved to: {summary_file}")


if __name__ == '__main__':
    main()
