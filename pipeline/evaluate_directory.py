"""
General evaluation script for running evaluation on predictions in a directory.

This script can evaluate predictions from any directory structure by:
1. Finding prediction files (TSV format)
2. Matching them with annotation files
3. Running genome_scores_evaluator on each pair
4. Aggregating and printing summary statistics

Usage:
    python evaluate_directory.py \
        --predictions_dir /path/to/predictions \
        --annotations_dir /path/to/annotations \
        --output_dir /path/to/output \
        --GAEF_dir /path/to/GAEF \
        [--threshold 0.5] \
        [--groovy_flag] \
        [--file_pattern "*.tsv"]
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
import os
import argparse
import json
import sys
import glob
import re
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'GAEF')))

from genome_scores_evaluator import main as genome_scores_evaluator, load_predictions_prop, load_annotations, evaluate_prediction_metrics
from GAEF.utils import Ontology


def verify_predictions_format(predictions_file: str, max_lines: int = 10) -> bool:
    """
    Verify predictions are in the correct format (optimized - only checks first few lines).
    
    Predictions file should be a TSV file with:
    - First column: protein_id
    - Remaining columns: GO:XXXXXX|score format
    
    Args:
        predictions_file: Path to predictions file
        max_lines: Maximum number of lines to check (default: 10)
        
    Returns:
        True if format is valid, False otherwise
    """
    try:
        with open(predictions_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line_num > max_lines:
                    break
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    return False
                for part in parts[1:]:
                    if '|' not in part:
                        return False
                    try:
                        float(part.split('|')[1])
                    except (ValueError, IndexError):
                        return False
        return True
    except Exception as _:
        return False


def find_prediction_files(predictions_dir: str, file_pattern: str = "*.tsv") -> List[str]:
    """
    Find all prediction files in the directory.
    
    Args:
        predictions_dir: Directory containing prediction files
        file_pattern: Glob pattern to match prediction files
        
    Returns:
        List of prediction file paths
    """
    pattern = os.path.join(predictions_dir, file_pattern)
    files = glob.glob(pattern)
    # Filter out non-prediction files if needed
    return sorted([f for f in files if os.path.isfile(f)])


def extract_identifier_from_filename(filename: str, pattern: Optional[str] = None) -> Optional[str]:
    """
    Extract identifier from filename for matching with annotation files.
    
    Supports multiple patterns:
    - taxon_<id> pattern
    - fold_<id>_taxon_<id> pattern
    - Custom regex pattern
    
    Args:
        filename: Name of the file (without path)
        pattern: Optional regex pattern to extract identifier
        
    Returns:
        Extracted identifier or None
    """
    if pattern:
        match = re.search(pattern, filename)
        if match:
            return match.group(1)
    
    # Try common patterns
    # Pattern 1: fold_<id>_taxon_<id>
    match = re.search(r'fold_(\d+)_taxon_([^_\.]+)', filename)
    if match:
        # return fold and taxon id
        return (match.group(1), match.group(2))

    # Pattern 2: taxon_<id>
    match = re.search(r'taxon_([^_\.]+)', filename)
    if match:
        return match.group(1)
    
    # Pattern 3: Extract base name without extension
    base_name = os.path.splitext(filename)[0]
    return base_name


def build_annotation_file_map(annotations_dir: str) -> Dict[str, str]:
    """
    Build a mapping from identifiers to annotation file paths for faster lookup.
    
    Args:
        annotations_dir: Directory containing annotation files
        
    Returns:
        Dictionary mapping identifier -> annotation file path
    """
    annotation_map = {}
    if not os.path.exists(annotations_dir):
        return annotation_map
    
    # Build map from existing files
    for filename in os.listdir(annotations_dir):
        if filename.endswith(('.tsv', '.goa')):
            # Try to extract identifier from filename
            identifier = extract_identifier_from_filename(filename)
            if identifier:
                filepath = os.path.join(annotations_dir, filename)
                # Allow multiple patterns to map to same file
                annotation_map[identifier] = filepath
                # Also add variations
                if filename.startswith('annots_taxon_'):
                    base_id = filename.replace('annots_taxon_', '').replace('.tsv', '').replace('.goa', '')
                    annotation_map[base_id] = filepath
                elif filename.startswith('annots_'):
                    base_id = filename.replace('annots_', '').replace('.tsv', '').replace('.goa', '')
                    annotation_map[base_id] = filepath
    
    return annotation_map


def find_matching_annotation_file(prediction_file: str, annotations_dir: str, 
                                  annotation_map: Optional[Dict[str, str]] = None,
                                  identifier: Optional[str] = None) -> Optional[str]:
    """
    Find matching annotation file for a prediction file.
    
    Args:
        prediction_file: Path to prediction file
        annotations_dir: Directory containing annotation files
        annotation_map: Pre-built mapping of identifiers to annotation files (for performance)
        identifier: Optional identifier extracted from prediction filename
        
    Returns:
        Path to annotation file or None if not found
    """
    if identifier is None:
        identifier = extract_identifier_from_filename(os.path.basename(prediction_file))
    
    if identifier is None:
        return None
    
    # Use pre-built map if available
    if annotation_map and identifier in annotation_map:
        return annotation_map[identifier]
    
    # Fallback to file system lookup
    patterns = [
        f'annots_taxon_{identifier}.tsv',
        f'annots_{identifier}.tsv',
        f'{identifier}.tsv',
        f'{identifier}.goa',
    ]
    
    for pattern in patterns:
        annotation_file = os.path.join(annotations_dir, pattern)
        if os.path.exists(annotation_file):
            return annotation_file
    
    return None


def evaluate_single_file(args: Tuple) -> Tuple[bool, str, Optional[str]]:
    """
    Evaluate a single prediction file (for parallel processing).
    
    Args:
        args: Tuple of (prediction_file, identifier, annotation_file, output_file, 
                       GAEF_dir, threshold, groovy_flag, verify_format, subontology)
        
    Returns:
        Tuple of (success, identifier, error_message)
    """
    (prediction_file, fold, identifier, annotation_file, output_file, 
     GAEF_dir, go_file, threshold, groovy_flag, verify_format, subontology) = args
    
    try:
        # Verify format if requested
        if verify_format and not verify_predictions_format(prediction_file):
            return (False, fold, identifier, f"Invalid format: {prediction_file}")
        
        # Run evaluation
        genome_scores_evaluator(
            assembly_name=identifier,
            annotations_file=annotation_file,
            predictions_file=prediction_file,
            GAEF_dir=GAEF_dir,
            go_file=go_file,
            groovy_flag=groovy_flag,
            python_consistency=not groovy_flag,
            python_ic=not groovy_flag,
            threshold=threshold,
            output_file=output_file,
            subontology=subontology,
        )
        return (True, fold, identifier, None)
    except Exception as e:
        return (False, fold, identifier, str(e))


def evaluate_directory(
    predictions_dir: str,
    annotations_dir: str,
    output_dir: str,
    GAEF_dir: str,
    go_file: str,
    threshold: float = 0.5,
    groovy_flag: bool = False,
    file_pattern: str = "*.tsv",
    identifier_pattern: Optional[str] = None,
    threshold_map: Optional[Dict[str, float]] = None,
    num_workers: Optional[int] = None,
    skip_existing: bool = True,
    verify_format: bool = True,
    subontology: str = 'cc',
) -> Dict:
    """
    Evaluate all predictions in a directory (optimized with parallel processing).
    
    Args:
        predictions_dir: Directory containing prediction files
        annotations_dir: Directory containing annotation files
        output_dir: Directory to save evaluation results
        GAEF_dir: Path to GAEF directory
        go_file: Path to GO file
        threshold: Default threshold for predictions
        groovy_flag: Use Groovy scripts for IC calculation
        file_pattern: Glob pattern to match prediction files
        identifier_pattern: Optional regex pattern to extract identifier from filenames
        threshold_map: Optional dictionary mapping identifiers to thresholds
        num_workers: Number of parallel workers (default: CPU count)
        skip_existing: Skip files that already have valid output files
        verify_format: Whether to verify prediction file format (can be slow)
        subontology: Subontology to evaluate (default: cc, choices: cc, bp, mf)
    Returns:
        Dictionary with summary statistics
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all prediction files
    prediction_files = find_prediction_files(predictions_dir, file_pattern)
    
    if not prediction_files:
        print(f"[Error] No prediction files found in {predictions_dir} with pattern {file_pattern}")
        return {}
    
    print(f"Found {len(prediction_files)} prediction files")
    
    # Build annotation file map for faster lookups
    print("Building annotation file map...")
    annotation_map = build_annotation_file_map(annotations_dir)
    print(f"Found {len(annotation_map)} annotation files")
    
    # Prepare evaluation tasks
    tasks = []
    skipped = 0
    errors_prep = []
    
    for prediction_file in prediction_files:
        # Extract identifier
        identifier = extract_identifier_from_filename(
            os.path.basename(prediction_file), 
            identifier_pattern
        )
        if isinstance(identifier, tuple):
            fold, taxon_id = identifier
            identifier = taxon_id
        else:
            fold = None
        
        if identifier is None:
            identifier = os.path.splitext(os.path.basename(prediction_file))[0]
        
        # Create output filename
        output_filename = f"evaluation_{fold}_{identifier}.json"
        output_file = os.path.join(output_dir, output_filename)
        
        # Skip if output already exists and is valid
        if skip_existing and os.path.exists(output_file):
            try:
                with open(output_file, 'r') as f:
                    json.load(f)  # Verify it's valid JSON
                skipped += 1
                continue
            except Exception as _:
                pass  # File exists but is invalid, re-evaluate
        
        # Get threshold for this file
        file_threshold = threshold
        if threshold_map and fold is not None and fold in threshold_map:
            file_threshold = threshold_map[fold]
        
        # Find matching annotation file
        annotation_file = find_matching_annotation_file(
            prediction_file, 
            annotations_dir, 
            annotation_map,
            identifier
        )
        
        if annotation_file is None:
            errors_prep.append(f"No annotation file found for {prediction_file} (fold: {fold}, identifier: {identifier})")
            continue
        
        # Add to tasks
        tasks.append((
            prediction_file,
            fold,
            identifier,
            annotation_file,
            output_file,
            GAEF_dir,
            go_file,
            file_threshold,
            groovy_flag,
            verify_format,
            subontology,
        ))
    
    if skipped > 0:
        print(f"Skipped {skipped} files that already have valid outputs")
    
    if errors_prep:
        print(f"\n[Warning] {len(errors_prep)} files skipped due to preparation errors:")
        for error in errors_prep[:5]:  # Show first 5
            print(f"  - {error}")
        if len(errors_prep) > 5:
            print(f"  ... and {len(errors_prep) - 5} more")
    
    if not tasks:
        print("No additional files to evaluate, computing micro-averaged metrics")
        summary = aggregate_metrics(output_dir, subontology)
        print_summary(summary, subontology)

        # Compute micro-averaged metrics
        micro_metrics = compute_micro_averaged_metrics(
            predictions_dir=predictions_dir,
            annotations_dir=annotations_dir,
            go_file=go_file,
            subontology=subontology,
            file_pattern=file_pattern,
            identifier_pattern=identifier_pattern,
        )
        print_summary(summary, subontology, micro_metrics)

        # Add micro metrics to summary
        if micro_metrics is not None:
            summary['micro_averaged_metrics'] = {subontology: micro_metrics}

        return summary
    
    # Determine number of workers
    if num_workers is None:
        num_workers = min(len(tasks), mp.cpu_count())
    
    print(f"\nEvaluating {len(tasks)} files with {num_workers} parallel workers...")
    
    # Evaluate files in parallel
    successful_evaluations = 0
    failed_evaluations = 0
    error_messages = []
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        future_to_task = {
            executor.submit(evaluate_single_file, task): task 
            for task in tasks
        }
        
        # Process results with progress bar
        with tqdm(total=len(tasks), desc="Evaluating predictions") as pbar:
            for future in as_completed(future_to_task):
                success, fold, identifier, error_msg = future.result()
                if success:
                    successful_evaluations += 1
                else:
                    failed_evaluations += 1
                    error_messages.append((fold, identifier, error_msg))
                pbar.update(1)
    
    # Print error summary
    if error_messages:
        print(f"\n[Errors] {len(error_messages)} evaluations failed:")
        for fold, identifier, error_msg in error_messages[:10]:  # Show first 10
            print(f"  - Fold {fold}: {identifier}: {error_msg}")
        if len(error_messages) > 10:
            print(f"  ... and {len(error_messages) - 10} more errors")
    
    print(f"\n{'='*60}")
    print("Evaluation Summary")
    print(f"{'='*60}")
    print(f"Total prediction files: {len(prediction_files)}")
    print(f"Skipped (already exist): {skipped}")
    print(f"Successful evaluations: {successful_evaluations}")
    print(f"Failed evaluations: {failed_evaluations}")
    
    # Aggregate and print summary statistics (macro-average)
    summary = aggregate_metrics(output_dir, subontology)
    
    # Compute micro-averaged metrics (CAFA standard)
    micro_metrics = compute_micro_averaged_metrics(
        predictions_dir=predictions_dir,
        annotations_dir=annotations_dir,
        go_file=go_file,
        subontology=subontology,
        file_pattern=file_pattern,
        identifier_pattern=identifier_pattern,
    )
    
    # Print both macro and micro averages
    print_summary(summary, subontology, micro_metrics)
    
    # Add micro metrics to summary
    if micro_metrics is not None:
        summary['micro_averaged_metrics'] = {subontology: micro_metrics}
    
    return summary


def aggregate_metrics(output_dir: str, subontology: str = 'cc') -> Dict:
    """
    Aggregate metrics from all evaluation JSON files.
    
    Args:
        output_dir: Directory containing evaluation JSON files
        subontology: Subontology to aggregate metrics for (default: cc)
        
    Returns:
        Dictionary with aggregated metrics
    """
    GAEF_metric_keys = [
        'essential_percentage', 
        'complete_has_part_percentage', 
        'complex_coherence', 
        'metacyc_complete_percentage'
    ]
    GAEF_metric_values = {key: [] for key in GAEF_metric_keys}
    taxon_consistent_values = []
    
    metric_keys = ['fmax', 'smin', 'avg_auc', 'aupr', 'precision', 'recall', 'tp', 'fp', 'fn']
    
    # Only aggregate metrics for the specified subontology
    ontology_metrics = {
        subontology: {key: [] for key in metric_keys},
    }
    
    json_files = [f for f in os.listdir(output_dir) if f.endswith('.json')]
    
    for json_file in json_files:
        json_path = os.path.join(output_dir, json_file)
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
                
                # Aggregate GAEF metrics
                for key in GAEF_metric_keys:
                    if key in data:
                        GAEF_metric_values[key].append(data[key])
                
                # Aggregate taxon consistency
                if 'satisfiable' in data:
                    taxon_consistent_values.append(data['satisfiable'])
                
                # Aggregate prediction metrics
                if 'prediction_metrics' in data:
                    # Per-ontology metrics (only for the specified subontology)
                    if subontology in data['prediction_metrics']:
                        for key in metric_keys:
                            if key in data['prediction_metrics'][subontology]:
                                ontology_metrics[subontology][key].append(data['prediction_metrics'][subontology][key])
        except Exception as e:
            print(f"[Warning] Error reading {json_file}: {e}")
            continue
    
    return {
        'GAEF_metrics': GAEF_metric_values,
        'taxon_consistency': taxon_consistent_values,
        'ontology_metrics': ontology_metrics,
        'total_files': len(json_files),
    }


def compute_micro_averaged_metrics(
    predictions_dir: str,
    annotations_dir: str,
    go_file: str,
    subontology: str = 'cc',
    file_pattern: str = "*.tsv",
    identifier_pattern: Optional[str] = None,
) -> Optional[Dict]:
    """
    Compute micro-averaged metrics by pooling all proteins from all genomes.
    This is the CAFA standard evaluation method.
    
    Args:
        predictions_dir: Directory containing prediction files
        annotations_dir: Directory containing annotation files
        go_file: Path to GO ontology file
        subontology: Subontology to evaluate (default: cc)
        file_pattern: Glob pattern to match prediction files
        identifier_pattern: Optional regex pattern to extract identifier from filenames
        
    Returns:
        Dictionary with micro-averaged metrics or None if computation fails
    """
    print(f"\n{'='*60}")
    print("Computing Micro-Averaged Metrics (CAFA Standard)")
    print(f"{'='*60}")
    
    # Find all prediction files
    prediction_files = find_prediction_files(predictions_dir, file_pattern)
    if not prediction_files:
        print("[Error] No prediction files found")
        return None
    
    # Build annotation file map
    annotation_map = build_annotation_file_map(annotations_dir)
    
    # Load ontology once
    go = Ontology(go_file)
    
    # Pool all predictions and annotations
    all_predictions = {}
    all_annotations = {}
    skipped_files = 0
    
    print(f"Loading and pooling {len(prediction_files)} prediction files...")
    for prediction_file in tqdm(prediction_files, desc="Loading files"):
        # Extract identifier
        identifier = extract_identifier_from_filename(
            os.path.basename(prediction_file), 
            identifier_pattern
        )
        if isinstance(identifier, tuple):
            _, identifier = identifier
        
        if identifier is None:
            identifier = os.path.splitext(os.path.basename(prediction_file))[0]
        
        # Find matching annotation file
        annotation_file = find_matching_annotation_file(
            prediction_file, 
            annotations_dir, 
            annotation_map,
            identifier
        )
        
        if annotation_file is None:
            skipped_files += 1
            continue
        
        try:
            # Load annotations for this genome
            genome_annotations = load_annotations(annotation_file)
            
            # Load predictions for this genome (only for annotated proteins)
            annotated_protein_ids = set(genome_annotations.keys())
            genome_predictions = load_predictions_prop(
                prediction_file, go, subontology, annotated_proteins=annotated_protein_ids
            )
            
            # Add to pooled data with genome prefix to avoid ID collisions
            for prot_id, go_terms in genome_annotations.items():
                pooled_id = f"{identifier}:{prot_id}"
                all_annotations[pooled_id] = go_terms
            
            for prot_id, scores in genome_predictions.items():
                pooled_id = f"{identifier}:{prot_id}"
                all_predictions[pooled_id] = scores
                
        except Exception as e:
            print(f"[Warning] Error loading {identifier}: {e}")
            skipped_files += 1
            continue
    
    if skipped_files > 0:
        print(f"[Warning] Skipped {skipped_files} files due to errors")
    
    print(f"\nPooled {len(all_annotations)} proteins from {len(prediction_files) - skipped_files} genomes")
    print("Computing micro-averaged F-max on entire pool...")
    
    # Compute metrics on pooled data
    try:
        metrics = evaluate_prediction_metrics(go, subontology, all_annotations, all_predictions)
        
        if metrics is not None:
            print(f"\nMicro-Averaged Metrics ({subontology.upper()}):")
            print(f"  Fmax: {metrics['fmax']:.5f}")
            print(f"  Smin: {metrics['smin']:.5f}")
            print(f"  Threshold: {metrics['tmax']:.5f}")
            print(f"  AUC: {metrics['avg_auc']:.5f}")
            print(f"  AUPR: {metrics['aupr']:.5f}")
            print(f"  AVGIC: {metrics['avgic']:.5f}")
            print(f"  Precision: {metrics['precision']:.5f}")
            print(f"  Recall: {metrics['recall']:.5f}")
            print(f"  TP: {metrics.get('tp', 'N/A')}, FP: {metrics.get('fp', 'N/A')}, FN: {metrics.get('fn', 'N/A')}")
            return metrics
        else:
            print("[Warning] Could not compute micro-averaged metrics")
            return None
            
    except Exception as e:
        print(f"[Error] Error computing micro-averaged metrics: {e}")
        import traceback
        traceback.print_exc()
        return None


def print_summary(summary: Dict, subontology: str = 'cc', micro_metrics: Optional[Dict] = None):
    """
    Print summary statistics.
    
    Args:
        summary: Dictionary with aggregated metrics (macro-average)
        subontology: Subontology to print metrics for (default: cc)
        micro_metrics: Optional dictionary with micro-averaged metrics (CAFA standard)
    """
    print(f"\n{'='*60}")
    print(f"Metrics Summary ({subontology.upper()} Ontology)")
    print(f"{'='*60}")
    print(f"Total evaluated files: {summary['total_files']}\n")
    
    # GAEF metrics
    print("--- GAEF Metrics (Macro-Average) ---")
    GAEF_metric_keys = [
        'essential_percentage', 
        'complete_has_part_percentage', 
        'complex_coherence', 
        'metacyc_complete_percentage'
    ]
    for key in GAEF_metric_keys:
        values = summary['GAEF_metrics'].get(key, [])
        if values:
            print(f"Average {key}: {np.mean(values):.3f}% (n={len(values)})")
    
    # Taxon consistency
    if summary['taxon_consistency']:
        print(f"Average taxon consistency: {np.mean(summary['taxon_consistency']) * 100:.3f}% (n={len(summary['taxon_consistency'])})")
    
    # Per-ontology metrics (macro-average: average per-genome)
    metric_keys = ['fmax', 'smin', 'avg_auc', 'aupr', 'precision', 'recall', 'tp', 'fp', 'fn']
    print(f"\n--- {subontology.upper()} Ontology Metrics (Macro-Average: Per-Genome) ---")
    if subontology in summary['ontology_metrics']:
        for key in metric_keys:
            values = summary['ontology_metrics'][subontology].get(key, [])
            if values:
                mean_val = np.mean(values)
                if key in ('tp', 'fp', 'fn'):
                    print(f"Average {key}: {mean_val:.1f} (n={len(values)})")
                else:
                    print(f"Average {key}: {mean_val:.5f} (n={len(values)})")
    
    # Micro-averaged metrics (CAFA standard: pooled proteins)
    if micro_metrics is not None:
        print(f"\n--- {subontology.upper()} Ontology Metrics (Micro-Average: CAFA Standard) ---")
        for key in metric_keys:
            if key in micro_metrics:
                if key in ('tp', 'fp', 'fn'):
                    print(f"{key}: {micro_metrics[key]}")
                else:
                    print(f"{key}: {micro_metrics[key]:.5f}")


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate predictions in a directory',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python evaluate_directory.py \\
      --predictions_dir ./predictions \\
      --annotations_dir ./annotations \\
      --output_dir ./evaluations \\
      --GAEF_dir ../../GAEF
  
  # With custom threshold file
  python evaluate_directory.py \\
      --predictions_dir ./predictions \\
      --annotations_dir ./annotations \\
      --output_dir ./evaluations \\
      --GAEF_dir ../../GAEF \\
      --threshold_file fold_thresholds.json
  
  # With custom file pattern
  python evaluate_directory.py \\
      --predictions_dir ./predictions \\
      --annotations_dir ./annotations \\
      --output_dir ./evaluations \\
      --GAEF_dir ../../GAEF \\
      --file_pattern "predictions_*.tsv"
        """
    )
    
    parser.add_argument('--predictions_dir', required=True,
                       help='Directory containing prediction files (TSV format)')
    parser.add_argument('--annotations_dir', required=True,
                       help='Directory containing annotation files')
    parser.add_argument('--output_dir', required=True,
                       help='Directory to save evaluation results (JSON files)')
    parser.add_argument('--GAEF_dir', default='../../GAEF',
                       help='Path to GAEF directory')
    parser.add_argument('--go_file', default="data/go-basic.obo", help='Path to GO file')  # 2025-10 version
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Default threshold for GO term scores (default: 0.5)')
    parser.add_argument('--threshold_file', type=str, default=None,
                       help='JSON file mapping identifiers to thresholds (e.g., {"taxon_123": 0.3})')
    parser.add_argument('--groovy_flag', action='store_true',
                       help='Use Groovy scripts for IC calculation and taxonomic consistency')
    parser.add_argument('--file_pattern', type=str, default='*.tsv',
                       help='Glob pattern to match prediction files (default: *.tsv)')
    parser.add_argument('--identifier_pattern', type=str, default=None,
                       help='Regex pattern to extract identifier from filenames (e.g., r"taxon_([^_\\.]+)")')
    parser.add_argument('--num_workers', type=int, default=None,
                       help='Number of parallel workers (default: CPU count)')
    parser.add_argument('--skip_existing', action='store_true',
                       help='Skip files that already have valid output files')
    parser.add_argument('--verify_format', action='store_true',
                       help='Verify format of prediction files')
    parser.add_argument('--subontology', type=str, default='cc', choices=['cc', 'bp', 'mf'],
                       help='Subontology to evaluate (default: cc, choices: cc, bp, mf)')
    
    args = parser.parse_args()
    
    # Load threshold map if provided
    threshold_map = None
    if args.threshold_file:
        if os.path.exists(args.threshold_file):
            with open(args.threshold_file, 'r') as f:
                threshold_map = json.load(f)
        else:
            print(f"[Warning] Threshold file {args.threshold_file} not found, using default threshold")
    
    # Run evaluation
    summary = evaluate_directory(
        predictions_dir=args.predictions_dir,
        annotations_dir=args.annotations_dir,
        output_dir=args.output_dir,
        GAEF_dir=args.GAEF_dir,
        go_file=args.go_file,
        threshold=args.threshold,
        groovy_flag=args.groovy_flag,
        file_pattern=args.file_pattern,
        identifier_pattern=args.identifier_pattern,
        threshold_map=threshold_map,
        num_workers=args.num_workers,
        skip_existing=args.skip_existing,
        verify_format=args.verify_format,
        subontology=args.subontology,
    )
    
    # Save summary to file
    summary_file = os.path.join(args.output_dir, 'evaluation_summary.json')
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary saved to: {summary_file}")


if __name__ == '__main__':
    main()

