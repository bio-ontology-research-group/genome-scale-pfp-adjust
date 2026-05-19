"""
This script processes proteome folds, holding one out a as a test set and the rest as a training set.
The annotation is done using Diamond-blast. 
A given list of model organisms are used as information available to all training sets but are not used for evaluation.
We have two sources of data:
- Swiss-Prot
- UniProt reference proteomes

For each held out test fold, we do the following:
- Run sequence similarity for all UniProt proteome proteins against the training set+model organisms
- Evaluate F_max of the annotations for taxa of test set


"""

import click as ck
import os
import glob
import json
import subprocess
import multiprocessing as mp
import gzip
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Set
from collections import defaultdict


import numpy as np
import pandas as pd
from tqdm import tqdm
import gc

try:
    import psutil
except ImportError:
    psutil = None

try:
    from scipy.sparse import csr_matrix, lil_matrix
    SCIPY_AVAILABLE = True
    print("[INFO] scipy.sparse is available")
except ImportError:
    SCIPY_AVAILABLE = False
    print("[WARNING] scipy.sparse not available. Memory usage will be significantly higher. "
          "Install scipy for optimal memory efficiency: pip install scipy")

def convert_predictions_to_tsv(predictions: Dict[str, Dict[str, float]], output_file: str, sparsity_threshold: float = 0.0) -> str:
    """
    Convert predictions dictionary to TSV format.

    Format per line:
      protein_id\tGO:XXXXXXX|score\tGO:YYYYYYY|score...
    Only include terms with score >= sparsity_threshold. Scores printed with 6 decimals.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        for protein_id, go_scores in predictions.items():
            if not go_scores:
                f.write(f"{protein_id}\n")
                continue
            parts = [protein_id]
            for go_term, score in go_scores.items():
                if score >= sparsity_threshold:
                    parts.append(f"{go_term}|{score:.6f}")
            f.write("\t".join(parts) + "\n")
    return output_file


def organize_predictions_by_taxon(predictions: Dict[str, Dict[str, float]],
                                  protein_to_taxid: Dict[str, str],
                                  allowed_taxa: Set[str]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Organize predictions by taxon, restricted to allowed_taxa.
    Returns: taxon_id -> {protein_id -> {GO_id: score}}
    """
    predictions_by_taxon: Dict[str, Dict[str, Dict[str, float]]] = {}
    missing_taxon = 0
    for protein_id, go_scores in predictions.items():
        taxon_id = protein_to_taxid.get(protein_id)
        if not taxon_id:
            missing_taxon += 1
            continue
        if allowed_taxa and taxon_id not in allowed_taxa:
            continue
        if taxon_id not in predictions_by_taxon:
            predictions_by_taxon[taxon_id] = {}
        predictions_by_taxon[taxon_id][protein_id] = go_scores
    ck.echo(f"[INFO] Organized predictions by taxon: {len(predictions_by_taxon)} taxa; missing taxon for {missing_taxon} proteins")
    return predictions_by_taxon


def write_taxon_tsvs(predictions_by_taxon: Dict[str, Dict[str, Dict[str, float]]],
                     base_dir: str,
                     split: str,
                     fold_id: int,
                     sparsity_threshold: float = 0.0) -> List[str]:
    """
    Write TSVs under {base_dir}/predictions/{split}/predictions_fold_{fold}_taxon_{taxid}.tsv
    """
    created_files: List[str] = []
    target_dir = os.path.join(base_dir, 'predictions', split)
    os.makedirs(target_dir, exist_ok=True)
    for taxon_id, proteins in predictions_by_taxon.items():
        if not proteins:
            continue
        output_file = os.path.join(target_dir, f"predictions_fold_{fold_id:02d}_taxon_{taxon_id}.tsv")
        convert_predictions_to_tsv(proteins, output_file, sparsity_threshold)
        created_files.append(output_file)
    ck.echo(f"[INFO] Wrote {len(created_files)} {split} prediction files to {target_dir}")
    return created_files

def load_model_organisms(model_organisms_file: str) -> Set[str]:
    """Load model organism taxids from file."""
    if not os.path.exists(model_organisms_file):
        raise FileNotFoundError(f"Model organisms file not found: {model_organisms_file}")
    
    model_organisms = set()
    try:
        with open(model_organisms_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                taxid = line.strip()
                if taxid and not taxid.startswith('#'):  # Skip comments
                    if not taxid.isdigit():
                        ck.echo(f"[WARNING] Invalid taxid '{taxid}' on line {line_num} in {model_organisms_file}")
                        continue
                    model_organisms.add(taxid)
    except IOError as e:
        raise IOError(f"Error reading model organisms file {model_organisms_file}: {e}")
    
    # if not model_organisms:
    #     raise ValueError(f"No valid model organism taxids found in {model_organisms_file}")
    print(f"[INFO] Loaded {len(model_organisms)} model organisms from {model_organisms_file}")
    
    return model_organisms


def load_proteome_mapping(proteomes_ids_file: str) -> Dict[str, str]:
    """
    Load taxon_id -> proteome_id mapping from uniprot_proteomes_ids.tsv file.
    
    Args:
        proteomes_ids_file: Path to uniprot_proteomes_ids.tsv file
        
    Returns:
        Dictionary mapping taxon_id to proteome_id
    """
    if not os.path.exists(proteomes_ids_file):
        raise FileNotFoundError(f"Proteomes IDs file not found: {proteomes_ids_file}")
    
    proteome_mapping = {}
    try:
        with open(proteomes_ids_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('TaxonID'):
                    continue  # Skip header
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    taxon_id = parts[0].strip()
                    proteome_id = parts[1].strip()
                    if taxon_id and proteome_id:
                        proteome_mapping[taxon_id] = proteome_id
    except IOError as e:
        raise IOError(f"Error reading proteomes IDs file {proteomes_ids_file}: {e}")
    
    ck.echo(f"[INFO] Loaded {len(proteome_mapping)} proteome mappings from {proteomes_ids_file}")
    return proteome_mapping


def load_annotations_from_tsv(annotations_dir: str, taxon_id: str) -> Dict[str, List[str]]:
    """
    Load annotations from TSV file for a given taxon.
    
    Args:
        annotations_dir: Directory containing annotation TSV files
        taxon_id: Taxon ID
        
    Returns:
        Dictionary mapping protein_id to list of GO terms
    """
    annotation_file = os.path.join(annotations_dir, f'annots_taxon_{taxon_id}.tsv')
    if not os.path.exists(annotation_file):
        return {}
    
    annotations = {}
    try:
        with open(annotation_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    protein_id = parts[0].strip()
                    go_terms = [term.strip() for term in parts[1:]] if len(parts) > 1 else []
                    if protein_id:
                        annotations[protein_id] = go_terms
    except IOError as e:
        ck.echo(f"[WARNING] Error reading annotation file {annotation_file}: {e}")
        return {}
    
    return annotations


def load_sequences_from_fasta(main_proteomes_dir: str, proteome_mapping: Dict[str, str], 
                               taxon_id: str) -> Dict[str, str]:
    """
    Load sequences from FASTA file for a given taxon.
    
    Args:
        main_proteomes_dir: Base directory for UniProt reference proteomes
        proteome_mapping: Dictionary mapping taxon_id to proteome_id
        taxon_id: Taxon ID
        
    Returns:
        Dictionary mapping protein_id to sequence
    """
    proteome_id = proteome_mapping.get(taxon_id)
    if not proteome_id:
        return {}
    
    proteome_file = find_proteome_fasta_file(main_proteomes_dir, proteome_id, taxon_id)
    if not proteome_file:
        return {}
    
    return extract_sequences_from_fasta(proteome_file)


def create_dataframe_from_annotations_and_sequences(annotations_dir: str, main_proteomes_dir: str,
                                                     proteome_mapping: Dict[str, str],
                                                     taxon_ids: Set[str]) -> pd.DataFrame:
    """
    Create a DataFrame from annotation TSV files and FASTA sequences for given taxon IDs.
    
    Args:
        annotations_dir: Directory containing annotation TSV files
        main_proteomes_dir: Base directory for UniProt reference proteomes
        proteome_mapping: Dictionary mapping taxon_id to proteome_id
        taxon_ids: Set of taxon IDs to process
        
    Returns:
        DataFrame with columns: proteins, sequences, prop_annotations, orgs
    """
    all_proteins = []
    all_sequences = []
    all_annotations = []
    all_orgs = []
    
    for taxon_id in taxon_ids:
        # Load annotations
        annotations = load_annotations_from_tsv(annotations_dir, taxon_id)
        
        # Load sequences
        sequences = load_sequences_from_fasta(main_proteomes_dir, proteome_mapping, taxon_id)
        
        # Combine annotations and sequences
        all_protein_ids = set(annotations.keys()) | set(sequences.keys())
        
        for protein_id in all_protein_ids:
            all_proteins.append(protein_id)
            all_sequences.append(sequences.get(protein_id, ''))
            all_annotations.append(annotations.get(protein_id, []))
            all_orgs.append(taxon_id)
    
    df = pd.DataFrame({
        'proteins': all_proteins,
        'sequences': all_sequences,
        'prop_annotations': all_annotations,
        'orgs': all_orgs
    })
    
    return df


def find_proteome_fasta_file(main_proteomes_dir: str, proteome_id: str, taxon_id: str) -> str:
    """
    Find the .fasta file for a given proteome_id and taxon_id.
    
    Args:
        main_proteomes_dir: Base directory for UniProt reference proteomes
        proteome_id: UniProt proteome ID (e.g., 'UP000005640')
        taxon_id: Taxon ID
        
    Returns:
        Path to .fasta file or None if not found
    """
    domains = ['Bacteria', 'Eukaryota', 'Archaea', 'Viruses']
    
    for domain in domains:
        proteome_file = os.path.join(main_proteomes_dir, domain, proteome_id, 
                                     f'{proteome_id}_{taxon_id}.fasta.gz')
        if os.path.exists(proteome_file):
            return proteome_file
    
    return None


def extract_sequences_from_fasta(proteome_file: str) -> Dict[str, str]:
    """
    Extract all protein sequences from a .fasta.gz file.
    
    Args:
        proteome_file: Path to .fasta.gz proteome file
        
    Returns:
        Dictionary mapping protein_id to sequence
    """
    
    sequences = {}

    with gzip.open(proteome_file, 'rt') as f:
        prot_id = ''
        seq = ''
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if prot_id and seq:
                    sequences[prot_id] = seq
                prot_id = line[1:].split('|')[2].split()[0]
                seq = ''
            elif line:
                seq += line
        
        if prot_id and seq:
            sequences[prot_id] = seq
    
    return sequences


def create_diamond_database(train_df: pd.DataFrame, model_organisms_df: pd.DataFrame, 
                          output_dir: str, fold_id: int) -> str:
    """
    Create Diamond database from training set + model organisms using protein sequences.
    
    Args:
        train_df: Training set proteins
        model_organisms_df: Model organism proteins
        output_dir: Output directory for Diamond database
        fold_id: Fold identifier
        
    Returns:
        Path to created Diamond database
    """
    if train_df.empty:
        raise ValueError("[ERROR] Training dataframe is empty")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        raise OSError(f"[ERROR] Failed to create output directory {output_dir}: {e}")
    
    # Create temporary FASTA file with training sequences
    fasta_file = os.path.join(output_dir, f'train_fold_{fold_id:02d}.fasta')
    db_file = os.path.join(output_dir, f'train_fold_{fold_id:02d}.dmnd')
    
    sequences_written = 0
    try:
        with open(fasta_file, 'w', encoding='utf-8') as f:
            # Add training set sequences
            for _, row in train_df.iterrows():
                protein_id = str(row.proteins)
                
                # Get sequence directly from DataFrame
                if hasattr(row, 'sequences') and row.sequences:
                    sequence = str(row.sequences).strip()
                    if sequence and sequence != 'nan':  # Skip empty or NaN sequences
                        f.write(f'>{protein_id}\n{sequence}\n')
                        sequences_written += 1
            
            # Add model organism sequences
            for _, row in model_organisms_df.iterrows():
                protein_id = str(row.proteins)
                
                # Get sequence directly from DataFrame
                if hasattr(row, 'sequences') and row.sequences:
                    sequence = str(row.sequences).strip()
                    if sequence and sequence != 'nan':  # Skip empty or NaN sequences
                        f.write(f'>{protein_id}\n{sequence}\n')
                        sequences_written += 1
    except IOError as e:
        raise IOError(f"[ERROR] Failed to write FASTA file {fasta_file}: {e}")
    
    if sequences_written == 0:
        raise ValueError(f"[ERROR] No valid sequences found to create database for fold {fold_id}")
    
    ck.echo(f"[INFO] Wrote {sequences_written} sequences to {fasta_file}")
    
    # Create Diamond database with optimized settings
    # Find Diamond executable (check if it's in PATH or current directory)
    diamond_exe = 'diamond'
    if os.path.exists('./diamond'):
        diamond_exe = './diamond'
    elif not any(os.path.exists(os.path.join(path, 'diamond')) for path in os.environ.get('PATH', '').split(os.pathsep)):
        raise RuntimeError("[ERROR] Diamond executable not found. Please ensure Diamond is installed and accessible.")
    
    cmd = [
        diamond_exe, 'makedb', 
        '--in', fasta_file, 
        '-d', db_file.replace('.dmnd', ''),
        '--threads', str(min(mp.cpu_count(), 16))  # Use multiple threads for database creation
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"[ERROR] Failed to create Diamond database: {result.stderr}")
    
    return db_file


def run_diamond_blast(query_fasta: str, db_file: str, output_file: str, 
                     evalue: float = 1e-5, max_target_seqs: int = 100, 
                     threads: int = None, memory_limit: str = None, 
                     batch_size: int = 10000) -> bool:
    """
    Run Diamond BLAST search with optimized settings for large queries.
    
    Args:
        query_fasta: Query FASTA file
        db_file: Diamond database file
        output_file: Output file for results
        evalue: E-value threshold
        max_target_seqs: Maximum number of target sequences
        threads: Number of threads (default: auto-detect)
        memory_limit: Memory limit (e.g., '32G')
        batch_size: Number of sequences to process in each batch (for large queries)
        
    Returns:
        True if successful, False otherwise
    """
    if threads is None:
        threads = mp.cpu_count()  # Use all available CPUs
    
    # Check if query file is large enough to warrant batching
    try:
        import os
        file_size_mb = os.path.getsize(query_fasta) / (1024 * 1024)
        if file_size_mb > 500:  # If query file > 500MB, use batching
            ck.echo(f"[INFO] Large query file detected ({file_size_mb:.1f} MB). Using batching...")
            # Adjust batch_size based on file size for optimal performance
            # For very large files (>1GB), use smaller batches
            if file_size_mb > 1000:
                adjusted_batch_size = min(batch_size, 200000)
            else:
                adjusted_batch_size = batch_size

            batch_success = run_diamond_blast_batched(query_fasta, db_file, output_file,
                                                      evalue, max_target_seqs, threads,
                                                      memory_limit, adjusted_batch_size)
            if batch_success:
                return True
            else:
                ck.echo(f"[WARNING] Batching failed, falling back to regular Diamond BLAST (may take longer)...")
                # Note: This fallback may fail due to memory issues with very large files
                pass
    except Exception as e:
        ck.echo(f"[WARNING] Could not check file size for batching: {e}")
        pass  # Continue with regular processing if size check fails
    
    # Prepare tmp directory before building command
    tmp_dir = os.path.join(os.path.dirname(output_file), 'tmp')
    os.makedirs(tmp_dir, exist_ok=True)

    # Find Diamond executable (check if it's in PATH or current directory)
    diamond_exe = 'diamond'
    if os.path.exists('./diamond'):
        diamond_exe = './diamond'
    elif not any(os.path.exists(os.path.join(path, 'diamond')) for path in os.environ.get('PATH', '').split(os.pathsep)):
        raise RuntimeError("[ERROR] Diamond executable not found. Please ensure Diamond is installed and accessible.")
    
    cmd = [
        diamond_exe, 'blastp',
        '-q', query_fasta,
        '-d', db_file.replace('.dmnd', ''),
        '-o', output_file,
        '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue',
        '--evalue', str(evalue),
        '--max-target-seqs', str(max_target_seqs),
        '--threads', str(threads),
        '--sensitive',           # Changed from --more-sensitive for ~2-3x speedup
        '--block-size', '16.0',  # Increased block size for better performance
        '--index-chunks', '16',  # Increased index chunks for better parallelization
        # '--compress', '1',       # Compress output to reduce I/O - DISABLED to prevent .txt.gz files
        '-t', tmp_dir,           # Use dedicated tmp directory (supported flag)
        '--min-score', '50',     # Minimum score threshold to reduce computation
        '--gapopen', '11',       # Optimized gap penalties
        '--gapextend', '1',      # Optimized gap extension
        '--matrix', 'BLOSUM62',  # Use standard scoring matrix
    ]

    # Add memory limit if specified
    if memory_limit:
        # Diamond accepts memory limit in GB (e.g., '64' for 64GB)
        # Extract numeric value if format like '64G' is provided
        mem_value = memory_limit.rstrip('GgBb').strip()
        try:
            # Validate it's a number
            float(mem_value)
            # Note: --memory-limit may not be available in all Diamond versions
            # Only add if it's a valid number
            cmd.extend(['--memory-limit', mem_value])
        except ValueError:
            pass  # Skip if not a valid number
    
    # Validate query FASTA file before running Diamond
    if not os.path.exists(query_fasta):
        ck.echo(f"[ERROR] Query FASTA file not found: {query_fasta}")
        return False
    
    if os.path.getsize(query_fasta) == 0:
        ck.echo(f"[ERROR] Query FASTA file is empty: {query_fasta}")
        return False
    
    # Check that the file starts with a FASTA header
    with open(query_fasta, 'r') as f:
        first_line = f.readline().strip()
        if not first_line.startswith('>'):
            ck.echo(f"[ERROR] Query FASTA file doesn't start with FASTA header: {query_fasta}")
            ck.echo(f"First line: '{first_line}'")
            return False
    
    ck.echo(f"[INFO] Running Diamond BLAST with {threads} threads...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        ck.echo(f"[ERROR] Diamond BLAST failed: {result.stderr}")
        return False
    
    return True


def run_diamond_blast_batched(query_fasta: str, db_file: str, output_file: str,
                             evalue: float, max_target_seqs: int, threads: int,
                             memory_limit: str, batch_size: int) -> bool:
    """
    Run Diamond BLAST search with batching for large query files.
    Splits the query file into smaller batches and processes them sequentially.

    Args:
        query_fasta: Query FASTA file
        db_file: Diamond database file
        output_file: Output file for results
        evalue: E-value threshold
        max_target_seqs: Maximum number of target sequences
        threads: Number of threads
        memory_limit: Memory limit
        batch_size: Number of sequences per batch

    Returns:
        True if successful, False otherwise
    """
    import tempfile
    import shutil

    ck.echo(f"[INFO] Splitting query file into batches of {batch_size} sequences...")

    # Create temporary directory for batch files
    tmp_dir = os.path.join(os.path.dirname(output_file), 'tmp_batches')
    os.makedirs(tmp_dir, exist_ok=True)

    try:
        # Split query FASTA into batches
        batch_files = []
        current_batch = []
        current_count = 0
        batch_num = 0

        with open(query_fasta, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_count >= batch_size:
                        # Write current batch
                        batch_file = os.path.join(tmp_dir, f'batch_{batch_num:04d}.fasta')
                        with open(batch_file, 'w') as bf:
                            bf.writelines(current_batch)
                        batch_files.append(batch_file)
                        current_batch = []
                        current_count = 0
                        batch_num += 1
                    current_count += 1
                current_batch.append(line)

        # Write last batch
        if current_batch:
            batch_file = os.path.join(tmp_dir, f'batch_{batch_num:04d}.fasta')
            with open(batch_file, 'w') as bf:
                bf.writelines(current_batch)
            batch_files.append(batch_file)

        ck.echo(f"[INFO] Created {len(batch_files)} batches")

        # Process each batch
        diamond_exe = 'diamond'
        if os.path.exists('./diamond'):
            diamond_exe = './diamond'

        # Create main tmp directory for Diamond temp files
        diamond_tmp = os.path.join(os.path.dirname(output_file), 'tmp')
        os.makedirs(diamond_tmp, exist_ok=True)

        # Clear output file
        with open(output_file, 'w') as f:
            pass

        for i, batch_file in enumerate(batch_files):
            ck.echo(f"[INFO] Processing batch {i+1}/{len(batch_files)}...")

            batch_output = os.path.join(tmp_dir, f'output_{i:04d}.txt')

            cmd = [
                diamond_exe, 'blastp',
                '-q', batch_file,
                '-d', db_file.replace('.dmnd', ''),
                '-o', batch_output,
                '--outfmt', '6', 'qseqid', 'sseqid', 'bitscore', 'evalue',
                '--evalue', str(evalue),
                '--max-target-seqs', str(max_target_seqs),
                '--threads', str(threads),
                '--sensitive',
                '--block-size', '150.0', # increased from 12.0 to 150.0 for better performance
                '--index-chunks', '4', # decreased from 16 to 4 for better parallelization
                '-t', diamond_tmp,
                '--min-score', '50',
                '--gapopen', '11',
                '--gapextend', '1',
                '--matrix', 'BLOSUM62',
            ]

            if memory_limit:
                mem_value = memory_limit.rstrip('GgBb').strip()
                try:
                    float(mem_value)
                    cmd.extend(['--memory-limit', mem_value])
                except ValueError:
                    pass

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode != 0:
                ck.echo(f"[ERROR] Diamond BLAST failed for batch {i+1}: {result.stderr}")
                return False

            # Append batch results to main output file
            if os.path.exists(batch_output):
                with open(output_file, 'a') as out_f:
                    with open(batch_output, 'r') as batch_f:
                        out_f.write(batch_f.read())

        ck.echo("[INFO] All batches processed successfully")
        return True

    finally:
        # Clean up temporary files
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)


def process_taxon_fasta(taxon_info: Tuple[str, str, str, Dict[str, str]]) -> int:
    """
    Process a single taxon's proteome file and write sequences to output FASTA.
    
    Args:
        taxon_info: Tuple of (taxid, main_proteomes_dir, output_file, proteome_mapping)
        
    Returns:
        Number of sequences processed
    """
    taxid, main_proteomes_dir, output_file, proteome_mapping = taxon_info
    sequences_written = 0
    
    if not taxid or not main_proteomes_dir or not output_file:
        return 0
    
    # Get proteome_id from taxon_id mapping
    proteome_id = proteome_mapping.get(taxid)
    if not proteome_id:
        return 0
    
    # Find the proteome file
    proteome_file = find_proteome_fasta_file(main_proteomes_dir, proteome_id, taxid)
    if not proteome_file:
        return 0
    
    # Extract sequences from the proteome file
    sequences = extract_sequences_from_fasta(proteome_file)
    if not sequences:
        return 0
    
    # Use a lock file to prevent concurrent writes to the same output file
    temp_file = output_file + f'.tmp.{taxid}'
    
    try:
        with open(temp_file, 'w', encoding='utf-8') as f:
            for protein_id, sequence in sequences.items():
                if sequence:  # Only write non-empty sequences
                    f.write(f'>{protein_id}\n{sequence}\n')
                    sequences_written += 1
        
        # Append temp file content to main output file with proper locking
        if sequences_written > 0:
            try:
                with open(output_file, 'a', encoding='utf-8') as out_f:
                    with open(temp_file, 'r', encoding='utf-8') as temp_f:
                        out_f.write(temp_f.read())
            except Exception:
                sequences_written = 0
    
    except Exception:
        sequences_written = 0
    finally:
        # Clean up temp file
        try:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        except Exception:
            pass
    
    return sequences_written


def add_proteins_to_query_fasta(query_fasta: str, output_fasta: str, proteins_df: pd.DataFrame) -> None:
    """
    Add protein sequences from DataFrame to the query FASTA file.
    
    Args:
        query_fasta: Path to existing query FASTA file
        output_fasta: Path to output combined FASTA file
        proteins_df: DataFrame containing proteins with 'proteins' and 'sequences' columns
    """
    sequences_added = 0
    
    try:
        with open(output_fasta, 'w', encoding='utf-8') as out_f:
            # Copy existing query sequences
            if os.path.exists(query_fasta):
                with open(query_fasta, 'r', encoding='utf-8') as query_f:
                    out_f.write(query_f.read())
            
            # Add protein sequences from DataFrame
            for _, row in proteins_df.iterrows():
                protein_id = str(row.proteins)

                # Get sequence directly from DataFrame
                if hasattr(row, 'sequences') and row.sequences:
                    sequence = str(row.sequences).strip()
                    if sequence and sequence != 'nan':  # Skip empty or NaN sequences
                        out_f.write(f'>{protein_id}\n{sequence}\n')
                        sequences_added += 1
    
    except IOError as e:
        raise IOError(f"[ERROR] Failed to write combined FASTA file {output_fasta}: {e}")
    
    ck.echo(f"[INFO] Added {sequences_added} protein sequences to query FASTA")


def create_query_fasta(test_organisms: Set[str], main_proteomes_dir: str, 
                      output_dir: str, fold_id: int, proteome_mapping: Dict[str, str],
                      max_workers: int = None) -> str:
    """
    Create query FASTA file for test organisms using parallel processing.
    
    Args:
        test_organisms: Set of test organism taxids
        main_proteomes_dir: Base directory for UniProt reference proteomes
        output_dir: Output directory
        fold_id: Fold identifier
        proteome_mapping: Dictionary mapping taxon_id to proteome_id
        max_workers: Maximum number of parallel workers
        
    Returns:
        Path to query FASTA file
    """
    os.makedirs(output_dir, exist_ok=True)
    query_fasta = os.path.join(output_dir, f'query_fold_{fold_id:02d}.fasta')
    
    # Clear the output file
    with open(query_fasta, 'w') as f:
        pass
    
    if max_workers is None:
        max_workers = min(len(test_organisms), mp.cpu_count())
    
    # Prepare arguments for parallel processing
    taxon_args = []
    for taxid in test_organisms:
        temp_file = os.path.join(output_dir, f'temp_{taxid}.fasta')
        taxon_args.append((taxid, main_proteomes_dir, temp_file, proteome_mapping))
    
    ck.echo(f"[INFO] Processing {len(test_organisms)} organisms with {max_workers} workers...")
    
    # Process organisms in parallel
    total_sequences = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(tqdm(
            executor.map(process_taxon_fasta, taxon_args),
            total=len(taxon_args),
            desc="Processing organisms"
        ))
        total_sequences = sum(results)
    
    # Combine all temporary files
    ck.echo("[INFO] Combining FASTA files...")
    valid_taxa = 0
    with open(query_fasta, 'w') as out_f:
        for taxid in test_organisms:
            temp_file = os.path.join(output_dir, f'temp_{taxid}.fasta')
            if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
                # Validate the temp file has at least one FASTA entry
                with open(temp_file, 'r') as in_f:
                    content = in_f.read().strip()
                    if content and content.startswith('>'):
                        out_f.write(content + '\n')
                        valid_taxa += 1
                os.remove(temp_file)  # Clean up temp file
            else:
                if os.path.exists(temp_file):
                    os.remove(temp_file)  # Clean up empty temp file
    
    ck.echo(f"[INFO] Combined FASTA files from {valid_taxa} organisms")
    
    ck.echo(f"[INFO] Created query FASTA with {total_sequences} sequences")
    return query_fasta


def get_protein_to_taxid_mapping(main_proteomes_dir: str, organism_taxids: Set[str], 
                                 proteome_mapping: Dict[str, str]) -> Dict[str, str]:
    """
    Create a mapping from UniProt protein IDs to organism taxids by parsing .dat.gz files.
    
    Args:
        main_proteomes_dir: Base directory for UniProt reference proteomes
        organism_taxids: Set of organism taxids to process
        proteome_mapping: Dictionary mapping taxon_id to proteome_id
        
    Returns:
        Dictionary mapping protein_id to taxid
    """
    protein_to_taxid = {}
    
    for taxid in organism_taxids:
        # Get proteome_id from taxon_id mapping
        proteome_id = proteome_mapping.get(taxid)
        if not proteome_id:
            continue
        
        # Find the proteome file
        proteome_file = find_proteome_fasta_file(main_proteomes_dir, proteome_id, taxid)
        if not proteome_file:
            continue
        
        # Extract sequences (which includes protein IDs) from the proteome file
        sequences = extract_sequences_from_fasta(proteome_file)
        
        # Map each protein ID to the taxid
        for protein_id in sequences.keys():
            protein_to_taxid[protein_id] = taxid
    
    ck.echo(f"[INFO] Mapped {len(protein_to_taxid)} proteins to {len(set(protein_to_taxid.values()))} organisms")
    return protein_to_taxid


def parse_diamond_results_chunk(chunk_lines: List[str]) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """
    Parse a chunk of Diamond BLAST results with optimized processing.
    
    Args:
        chunk_lines: List of lines from Diamond output
        
    Returns:
        Dictionary mapping query protein to subject proteins and (bitscore, evalue) tuples
    """
    results = defaultdict(dict)
    
    # Pre-compile patterns for faster processing
    import re
    tab_pattern = re.compile(r'\t')
    
    for line in chunk_lines:
        if not line.strip():  # Skip empty lines
            continue
            
        parts = tab_pattern.split(line.strip(), 3)  # Split only first 3 tabs
        if len(parts) >= 4:
            query_id = parts[0]
            subject_id = parts[1]
            try:
                bitscore = float(parts[2])
                evalue = float(parts[3])
                
                # Store best hit per query-subject pair (based on bitscore)
                if subject_id not in results[query_id] or bitscore > results[query_id][subject_id][0]:
                    results[query_id][subject_id] = (bitscore, evalue)
            except (ValueError, IndexError):
                continue  # Skip malformed lines
    
    return results


def parse_diamond_results(diamond_file: str, chunk_size: int = 1000000) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """
    Parse Diamond BLAST results with optimized chunked processing for large files.
    
    Args:
        diamond_file: Path to Diamond output file
        chunk_size: Number of lines to process in each chunk
        
    Returns:
        Dictionary mapping query protein to subject proteins and (bitscore, evalue) tuples
    """
    results = defaultdict(dict)
    
    if not os.path.exists(diamond_file):
        ck.echo(f"[WARNING] Diamond results file not found: {diamond_file}")
        return results
    
    ck.echo(f"[INFO] Parsing Diamond results from {diamond_file}...")
    
    try:
        # Get file size for progress estimation
        file_size = os.path.getsize(diamond_file)
        ck.echo(f"[INFO] File size: {file_size / (1024**3):.2f} GB")
        
        # Adaptive chunk size based on available memory
        if psutil is not None:
            try:
                available_memory_gb = psutil.virtual_memory().available / (1024**3)
                if available_memory_gb < 4:  # Less than 4GB available
                    chunk_size = min(chunk_size, 100000)
                    ck.echo(f"[INFO] Reducing chunk size to {chunk_size} due to limited memory")
            except Exception:
                pass  # Continue with default chunk size if memory check fails
        
        # Use buffered reading for better I/O performance
        with open(diamond_file, 'r', buffering=65536, encoding='utf-8', errors='ignore') as f:
            chunk_lines = []
            processed_lines = 0
            
            for line in f:
                chunk_lines.append(line)
                
                if len(chunk_lines) >= chunk_size:
                    # Process chunk in parallel if large enough
                    if chunk_size >= 100000:  # Use parallel processing for large chunks
                        chunk_results = parse_diamond_results_chunk_parallel(chunk_lines)
                    else:
                        chunk_results = parse_diamond_results_chunk(chunk_lines)
                    
                    # Merge results efficiently
                    for query_id, hits in chunk_results.items():
                        if query_id not in results:
                            results[query_id] = hits
                        else:
                            for subject_id, (bitscore, evalue) in hits.items():
                                if subject_id not in results[query_id] or bitscore > results[query_id][subject_id][0]:
                                    results[query_id][subject_id] = (bitscore, evalue)
                    
                    processed_lines += len(chunk_lines)
                    if processed_lines % 1000000 == 0:  # Log every million lines
                        ck.echo(f"[INFO] Processed {processed_lines:,} lines...")
                        # Periodic memory cleanup
                        gc.collect()
                    
                    # Clear chunk to free memory
                    chunk_lines.clear()
                    del chunk_results  # Explicitly delete to free memory
            
            # Process remaining lines
            if chunk_lines:
                chunk_results = parse_diamond_results_chunk(chunk_lines)
                for query_id, hits in chunk_results.items():
                    if query_id not in results:
                        results[query_id] = hits
                    else:
                        for subject_id, (bitscore, evalue) in hits.items():
                            if subject_id not in results[query_id] or bitscore > results[query_id][subject_id][0]:
                                results[query_id][subject_id] = (bitscore, evalue)
    
    except (IOError, OSError) as e:
        ck.echo(f"[ERROR] Error reading Diamond results file {diamond_file}: {e}")
        return defaultdict(dict)
    except Exception as e:
        ck.echo(f"[ERROR] Unexpected error parsing Diamond results: {e}")
        return defaultdict(dict)
    
    ck.echo(f"[INFO] Parsed {len(results)} query proteins with hits")
    return results


def parse_diamond_results_chunk_parallel(chunk_lines: List[str]) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """
    Parse a chunk of Diamond BLAST results using parallel processing.
    
    Args:
        chunk_lines: List of lines from Diamond output
        
    Returns:
        Dictionary mapping query protein to subject proteins and (bitscore, evalue) tuples
    """
    from concurrent.futures import ThreadPoolExecutor
    import math
    
    # Split chunk into smaller sub-chunks for parallel processing
    num_workers = min(mp.cpu_count(), 8)  # Limit to 8 workers
    sub_chunk_size = max(1000, len(chunk_lines) // num_workers)
    
    sub_chunks = []
    for i in range(0, len(chunk_lines), sub_chunk_size):
        sub_chunks.append(chunk_lines[i:i + sub_chunk_size])
    
    results = defaultdict(dict)
    
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(parse_diamond_results_chunk, sub_chunk) for sub_chunk in sub_chunks]
        
        for future in futures:
            chunk_results = future.result()
            for query_id, hits in chunk_results.items():
                if query_id not in results:
                    results[query_id] = hits
                else:
                    for subject_id, (bitscore, evalue) in hits.items():
                        if subject_id not in results[query_id] or bitscore > results[query_id][subject_id][0]:
                            results[query_id][subject_id] = (bitscore, evalue)
    
    return results


def compute_diamond_predictions(diamond_results: Dict[str, Dict[str, Tuple[float, float]]], 
                               train_annotations: Dict[str, Set[str]], 
                               terms: List[str]) -> Dict[str, Dict[str, float]]:
    """
    Compute GO term predictions from Diamond results.
    
    Args:
        diamond_results: Diamond BLAST results with (bitscore, evalue) tuples
        train_annotations: Training protein annotations
        terms: List of GO terms
        
    Returns:
        Dictionary mapping protein to GO term predictions
    """
    predictions = {}
    
    if not diamond_results or not train_annotations or not terms or len(terms) == 0:
        ck.echo("[WARNING] Empty input data for prediction computation")
        return predictions
    
    term_index = {term: i for i, term in enumerate(terms)}
    n_terms = len(terms)
    
    # Process in batches to manage memory usage (larger batches for better performance)
    batch_size = 5000  # Increased from 1000 for better performance
    query_proteins = list(diamond_results.keys())
    
    for batch_start in tqdm(range(0, len(query_proteins), batch_size), desc="Computing predictions"):
        batch_end = min(batch_start + batch_size, len(query_proteins))
        batch_proteins = query_proteins[batch_start:batch_end]
        
        for query_protein in batch_proteins:
            hits = diamond_results[query_protein]
            if not hits:
                continue
            
            # Get all GO terms from hits and compute scores in one pass
            term_score_sums = defaultdict(float)
            total_score = 0.0
            
            for subject_protein, (bitscore, evalue) in hits.items():
                if subject_protein in train_annotations:
                    # Accumulate scores for all terms annotated to this subject
                    for term in train_annotations[subject_protein]:
                        if term in term_index:
                            term_score_sums[term] += bitscore
                    total_score += bitscore
            
            if total_score == 0 or not term_score_sums:
                continue
            
            # Normalize scores
            term_scores = {term: score_sum / total_score 
                          for term, score_sum in term_score_sums.items() 
                          if score_sum > 0}
            
            if term_scores:
                predictions[query_protein] = term_scores
        
        # Periodic memory cleanup for very large datasets
        if batch_start > 0 and batch_start % 50000 == 0:
            gc.collect()
    
    ck.echo(f"[INFO] Generated predictions for {len(predictions)} proteins")
    return predictions


