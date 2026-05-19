"""
Evaluate genome annotations against GAEF metrics (completeness, coherence, consistency, IC).

Unlike evaluate_directory.py which evaluates predictions (with scores) against annotations,
this script operates on raw annotation files directly — no prediction quality metrics
(fmax, smin, etc.) are computed.

Input format: TSV files with protein_id\tGO_term1\tGO_term2\t...

Usage:
    python evaluate_annotations_gaef.py \
        --annotations_dir /path/to/annotations \
        --output_dir /path/to/output \
        --GAEF_dir /path/to/GAEF \
        [--go_file data/go-basic.obo] \
        [--file_pattern "annots_taxon_*.tsv"] \
        [--num_workers 16] \
        [--skip_existing]
"""

import os
import sys
import argparse
import json
import glob
import re
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'GAEF')))

from genome_scores_evaluator import load_annotations, GAEF_evaluation_from_terms


def find_annotation_files(annotations_dir, file_pattern="annots_taxon_*.tsv"):
    """Find all annotation files matching the pattern in the directory."""
    pattern = os.path.join(annotations_dir, file_pattern)
    files = glob.glob(pattern)
    return sorted([f for f in files if os.path.isfile(f)])


def extract_taxon_id(filename):
    """
    Extract taxon ID from annotation filename.

    Supports patterns:
    - annots_taxon_<id>.tsv
    - annots_<id>.tsv
    - Falls back to base filename without extension
    """
    # Pattern 1: annots_taxon_<id>.tsv
    match = re.search(r'annots_taxon_([^_\.]+)', filename)
    if match:
        return match.group(1)

    # Pattern 2: annots_<id>.tsv
    match = re.search(r'annots_([^_\.]+)', filename)
    if match:
        return match.group(1)

    # Fallback: base name without extension
    return os.path.splitext(filename)[0]


def evaluate_single_annotation(args):
    """
    Worker function: evaluate GAEF metrics for a single annotation file.

    Args:
        args: Tuple of (annotation_file, taxon_id, output_file, GAEF_dir, go_file, subontology)

    Returns:
        Tuple of (success, taxon_id, error_message)
    """
    annotation_file, taxon_id, output_file, GAEF_dir, go_file, subontology = args

    try:
        protein_go_terms = load_annotations(annotation_file)
        if not protein_go_terms:
            return (False, taxon_id, f"No annotations loaded from {annotation_file}")

        GAEF_evaluation_from_terms(
            assembly_name=taxon_id,
            protein_go_terms=protein_go_terms,
            GAEF_dir=GAEF_dir,
            go_file=go_file,
            subontology=subontology,
            output_file=output_file,
        )
        return (True, taxon_id, None)
    except Exception as e:
        return (False, taxon_id, str(e))


def evaluate_annotations_directory(
    annotations_dir,
    output_dir,
    GAEF_dir,
    go_file,
    file_pattern="annots_taxon_*.tsv",
    num_workers=None,
    skip_existing=False,
    subontology=None,
):
    """
    Evaluate GAEF metrics for all annotation files in a directory.

    Args:
        annotations_dir: Directory containing annotation TSV files
        output_dir: Directory to save evaluation JSON files
        GAEF_dir: Path to GAEF directory
        go_file: Path to GO OBO file
        file_pattern: Glob pattern for annotation files
        num_workers: Number of parallel workers (default: CPU count)
        skip_existing: Skip files that already have valid JSON output
        subontology: Optional filter (cc/bp/mf, default: all)

    Returns:
        Dictionary with aggregated summary
    """
    os.makedirs(output_dir, exist_ok=True)

    # Find annotation files
    annotation_files = find_annotation_files(annotations_dir, file_pattern)
    if not annotation_files:
        print(f"[Error] No annotation files found in {annotations_dir} with pattern {file_pattern}")
        return {}

    print(f"Found {len(annotation_files)} annotation files")

    # Prepare tasks
    tasks = []
    skipped = 0

    for annotation_file in annotation_files:
        taxon_id = extract_taxon_id(os.path.basename(annotation_file))

        output_file = os.path.join(output_dir, f"gaef_{taxon_id}_report.json")

        # Skip if output already exists and is valid
        if skip_existing and os.path.exists(output_file):
            try:
                with open(output_file, 'r') as f:
                    json.load(f)
                skipped += 1
                continue
            except (json.JSONDecodeError, IOError):
                pass  # Invalid file, re-evaluate

        tasks.append((
            annotation_file,
            taxon_id,
            output_file,
            GAEF_dir,
            go_file,
            subontology,
        ))

    if skipped > 0:
        print(f"Skipped {skipped} files with existing valid outputs")

    if not tasks:
        print("No additional files to evaluate")
        summary = aggregate_gaef_metrics(output_dir)
        print_gaef_summary(summary, output_dir)
        return summary

    # Determine workers
    if num_workers is None:
        num_workers = min(len(tasks), mp.cpu_count())

    print(f"\nEvaluating {len(tasks)} files with {num_workers} parallel workers...")

    # Evaluate in parallel
    successful = 0
    failed = 0
    error_messages = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_task = {
            executor.submit(evaluate_single_annotation, task): task
            for task in tasks
        }

        with tqdm(total=len(tasks), desc="Evaluating annotations") as pbar:
            for future in as_completed(future_to_task):
                success, taxon_id, error_msg = future.result()
                if success:
                    successful += 1
                else:
                    failed += 1
                    error_messages.append((taxon_id, error_msg))
                pbar.update(1)

    # Print error summary
    if error_messages:
        print(f"\n[Errors] {len(error_messages)} evaluations failed:")
        for taxon_id, error_msg in error_messages[:10]:
            print(f"  - {taxon_id}: {error_msg}")
        if len(error_messages) > 10:
            print(f"  ... and {len(error_messages) - 10} more errors")

    print(f"\n{'='*60}")
    print(f"Evaluation Summary")
    print(f"{'='*60}")
    print(f"Total annotation files: {len(annotation_files)}")
    print(f"Skipped (already exist): {skipped}")
    print(f"Successful evaluations: {successful}")
    print(f"Failed evaluations: {failed}")

    # Aggregate and print
    summary = aggregate_gaef_metrics(output_dir)
    print_gaef_summary(summary, output_dir)

    return summary


def aggregate_gaef_metrics(output_dir):
    """
    Aggregate GAEF metrics from all JSON report files in the output directory.

    Returns:
        Dictionary with per-taxon records, aggregate lists per metric, and total file count.
    """
    metric_keys = [
        'essential_percentage',
        'complete_has_part_percentage',
        'complex_coherence',
        'metacyc_complete_percentage',
        'ic_depth',
        'ic_breadth',
        'normalized_ic_breadth',
    ]
    metrics = {key: [] for key in metric_keys}
    taxon_consistency = []
    per_taxon = []

    json_files = sorted([f for f in os.listdir(output_dir) if f.endswith('.json') and f != 'evaluation_summary.json'])

    for json_file in json_files:
        json_path = os.path.join(output_dir, json_file)
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)

            taxon_id = data.get('assembly_name', os.path.splitext(json_file)[0])
            record = {'taxon': taxon_id}

            for key in metric_keys:
                val = data.get(key)
                record[key] = val
                if val is not None:
                    metrics[key].append(val)

            satisfiable = data.get('satisfiable')
            record['satisfiable'] = satisfiable
            if satisfiable is not None:
                taxon_consistency.append(satisfiable)

            per_taxon.append(record)
        except Exception as e:
            print(f"[Warning] Error reading {json_file}: {e}")

    return {
        'per_taxon': per_taxon,
        'metrics': metrics,
        'taxon_consistency': taxon_consistency,
        'total_files': len(json_files),
    }


def print_gaef_summary(summary, output_dir=None):
    """Print per-taxon metrics table and aggregate summary with mean/std/min/max."""
    metric_keys = [
        'essential_percentage',
        'complete_has_part_percentage',
        'complex_coherence',
        'metacyc_complete_percentage',
        'ic_depth',
        'ic_breadth',
        'normalized_ic_breadth',
    ]
    short_labels = {
        'essential_percentage': 'Essential%',
        'complete_has_part_percentage': 'HasPart%',
        'complex_coherence': 'Complex%',
        'metacyc_complete_percentage': 'MetaCyc%',
        'ic_depth': 'IC_depth',
        'ic_breadth': 'IC_breadth',
        'normalized_ic_breadth': 'NormIC',
    }
    label_map = {
        'essential_percentage': 'Essential terms %',
        'complete_has_part_percentage': 'Process coherence (has-part %)',
        'complex_coherence': 'Complex coherence %',
        'metacyc_complete_percentage': 'MetaCyc pathway completion %',
        'ic_depth': 'IC depth',
        'ic_breadth': 'IC breadth',
        'normalized_ic_breadth': 'Normalized IC breadth',
    }

    per_taxon = summary.get('per_taxon', [])

    # Determine which metric columns have data
    active_keys = [k for k in metric_keys if any(r.get(k) is not None for r in per_taxon)]

    # Print per-taxon table
    if per_taxon:
        print(f"\n{'='*60}")
        print(f"Per-Taxon GAEF Metrics ({len(per_taxon)} taxa)")
        print(f"{'='*60}")

        # Header
        header_parts = [f"{'Taxon':<20s}", f"{'Consistent':>10s}"]
        header_parts += [f"{short_labels[k]:>12s}" for k in active_keys]
        header = "  ".join(header_parts)
        print(header)
        print("-" * len(header))

        # Rows sorted by taxon
        for record in sorted(per_taxon, key=lambda r: r['taxon']):
            parts = [f"{record['taxon']:<20s}"]
            sat = record.get('satisfiable')
            parts.append(f"{'Yes' if sat else 'No':>10s}" if sat is not None else f"{'N/A':>10s}")
            for k in active_keys:
                val = record.get(k)
                parts.append(f"{val:>12.4f}" if val is not None else f"{'N/A':>12s}")
            print("  ".join(parts))

        # Save per-taxon TSV
        if output_dir:
            tsv_path = os.path.join(output_dir, 'per_taxon_gaef_metrics.tsv')
            with open(tsv_path, 'w') as f:
                tsv_header = ['taxon', 'satisfiable'] + active_keys
                f.write('\t'.join(tsv_header) + '\n')
                for record in sorted(per_taxon, key=lambda r: r['taxon']):
                    row = [record['taxon']]
                    sat = record.get('satisfiable')
                    row.append(str(sat) if sat is not None else '')
                    for k in active_keys:
                        val = record.get(k)
                        row.append(f"{val:.6f}" if val is not None else '')
                    f.write('\t'.join(row) + '\n')
            print(f"\nPer-taxon TSV saved to: {tsv_path}")

    # Aggregate summary
    print(f"\n{'='*60}")
    print(f"Aggregate GAEF Metrics Summary")
    print(f"{'='*60}")
    print(f"Total evaluated files: {summary['total_files']}\n")

    for key in active_keys:
        values = summary['metrics'].get(key, [])
        if values:
            arr = np.array(values)
            print(f"{label_map[key]}:")
            print(f"  mean={arr.mean():.4f}  std={arr.std():.4f}  min={arr.min():.4f}  max={arr.max():.4f}  n={len(arr)}")

    # Taxon consistency
    tc = summary['taxon_consistency']
    if tc:
        consistent = sum(1 for v in tc if v)
        total = len(tc)
        print(f"\nTaxonomic consistency: {consistent}/{total} satisfiable ({consistent/total*100:.1f}%)")
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate genome annotations against GAEF metrics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python evaluate_annotations_gaef.py \\
      --annotations_dir /path/to/annotations \\
      --output_dir /path/to/output \\
      --GAEF_dir ../../GAEF

  python evaluate_annotations_gaef.py \\
      --annotations_dir /path/to/annotations \\
      --output_dir /path/to/output \\
      --GAEF_dir ../../GAEF \\
      --skip_existing \\
      --num_workers 16
        """
    )

    parser.add_argument('--annotations_dir', required=True,
                        help='Directory containing annotation TSV files')
    parser.add_argument('--output_dir', required=True,
                        help='Directory to save evaluation JSON reports')
    parser.add_argument('--GAEF_dir', default='../../GAEF',
                        help='Path to GAEF directory (default: ../../GAEF)')
    parser.add_argument('--go_file', default='data/go-basic.obo',
                        help='Path to GO OBO file (default: data/go-basic.obo)')
    parser.add_argument('--file_pattern', default='annots_taxon_*.tsv',
                        help='Glob pattern for annotation files (default: annots_taxon_*.tsv)')
    parser.add_argument('--num_workers', type=int, default=None,
                        help='Number of parallel workers (default: CPU count)')
    parser.add_argument('--skip_existing', action='store_true',
                        help='Skip files that already have valid output JSON')
    parser.add_argument('--subontology', type=str, default=None, choices=['cc', 'bp', 'mf'],
                        help='Subontology to evaluate (default: all)')

    args = parser.parse_args()

    summary = evaluate_annotations_directory(
        annotations_dir=args.annotations_dir,
        output_dir=args.output_dir,
        GAEF_dir=args.GAEF_dir,
        go_file=args.go_file,
        file_pattern=args.file_pattern,
        num_workers=args.num_workers,
        skip_existing=args.skip_existing,
        subontology=args.subontology,
    )

    # Save summary
    summary_file = os.path.join(args.output_dir, 'evaluation_summary.json')
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"Summary saved to: {summary_file}")


if __name__ == '__main__':
    main()
