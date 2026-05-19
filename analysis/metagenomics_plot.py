"""
# Plots motivated by metagenomics

The idea is to sample contigs of varying sizes (1000, 10000, 100000, 1000000) and obtain their predicted taxon assignment. Can also compare against Kraken2.
    - need to modify current adjust script to extract the predicted taxon assignment
"""

LENGTHS = [1000, 10000, 100000, 1000000]
NUM_SAMPLES = 50
SUBONTOLOGIES = ['cc', 'mf', 'bp']

import os
import sys
import csv
import argparse
import random
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAXON_CONSISTENCY_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'taxon_consistency'
)
sys.path.append(TAXON_CONSISTENCY_DIR)

from adjust_ortools import (
    load_predictions,
    load_constraints,
    load_go_hierarchy,
    load_taxon_hierarchy,
    compute_demotion_flip_cost,
    get_annotated_predictions,
    solve_sat_ortools_taxon_assignment,
    get_all_ancestors,
    normalize_taxon_id,
)

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DEFAULT_CONSTRAINTS = os.path.join(BASE_DIR, 'data', 'go_taxon_constraints_extracted_obo.tsv')
DEFAULT_GO_HIERARCHY = os.path.join(BASE_DIR, 'data', 'go_hierarchy.tsv')
DEFAULT_TAXON_HIERARCHY = os.path.join(BASE_DIR, 'data', 'taxon_hierarchy.tsv')
DEFAULT_NCBITAXON_HIERARCHY = os.path.join(BASE_DIR, 'data', 'ncbitaxon_hierarchy.tsv')


def is_matching_taxon(predicted_taxons_list, true_taxon_lineage):
    return len(set(predicted_taxons_list) - set(true_taxon_lineage)) == 0


def is_general_matching(is_matching, most_specific_taxon) -> bool:
    """A sample is 'generally matching' if:
    1. is_matching is True and no specific taxon was predicted (vacuous match), OR
    2. the most specific predicted taxon is NCBITaxon_131567 (cellular organisms).
    """
    if is_matching and not most_specific_taxon:
        return True
    if is_matching and most_specific_taxon == 'NCBITaxon_131567':
        return True
    return False


def find_predictions_file(predictions_dir: str, subontology: str, taxon_id: str) -> Optional[str]:
    """
    Search for a taxon's prediction TSV file under predictions_dir/<subontology>/.
    Returns the path of the first .tsv file whose name contains taxon_id, or None.
    """
    sub_dir = os.path.join(predictions_dir, subontology, "predictions")
    if not os.path.isdir(sub_dir):
        return None
    for fname in os.listdir(sub_dir):
        if taxon_id in fname and fname.endswith('.tsv'):
            return os.path.join(sub_dir, fname)
    return None


def load_cds_positions(cds_file: str) -> List[Tuple[int, int, str]]:
    """
    Load CDS positions from a TSV produced by metagenomics_find_locations.py.
    Column 4 (index 3) is expected to be the UniProt entry name directly.
    Returns list of (start, end, entry_name).
    """
    entries = []
    with open(cds_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader, None)  # skip header
        for row in reader:
            if len(row) < 4:
                continue
            start = int(row[0])
            end = int(row[1])
            entry_name = row[3]
            entries.append((start, end, entry_name))
    print(f"[INFO] Loaded {len(entries)} CDS entries from {cds_file}")
    return entries


def sample_contig_proteins(
    cds_entries: List[Tuple[int, int, str]],
    genome_size: int,
    contig_length: int,
    rng: random.Random,
) -> Tuple[int, List[str]]:
    """
    Sample a random contig and return proteins overlapping with it.
    Returns (contig_start, list_of_entry_names).
    """
    max_start = genome_size - contig_length
    if max_start < 1:
        return 1, [e[2] for e in cds_entries]
    contig_start = rng.randint(1, max_start)
    contig_end = contig_start + contig_length

    seen: set = set()
    overlapping = []
    for start, end, entry_name in cds_entries:
        if start <= contig_end and end >= contig_start and entry_name not in seen:
            seen.add(entry_name)
            overlapping.append(entry_name)
    return contig_start, overlapping


def get_most_specific_taxon(
    predicted_taxons: List[str],
    taxon_is_a_hierarchy: Dict[str, Set[str]],
) -> Optional[str]:
    """
    Return the most specific (deepest) predicted taxon, excluding Union terms.
    Depth = length of ancestor path to root.
    """
    ancestor_cache = {}
    best_taxon = None
    best_depth = -1
    for taxon in predicted_taxons:
        if taxon.startswith('NCBITaxon_Union_'):
            continue
        ancestors = get_all_ancestors(taxon, taxon_is_a_hierarchy, ancestor_cache)
        depth = len(ancestors)
        if depth > best_depth:
            best_depth = depth
            best_taxon = taxon
    return best_taxon


def run_metagenomics_experiment(
    predictions: Dict[str, Dict[str, float]],
    cds_entries: List[Tuple[int, int, str]],
    genome_size: int,
    in_taxon_constraints: Dict[str, List[str]],
    never_in_taxon_constraints: Dict[str, List[str]],
    go_hierarchy: Dict[str, Set[str]],
    taxon_is_a_hierarchy: Dict[str, Set[str]],
    taxon_disjoint_from_hierarchy: Dict[str, Set[str]],
    taxon_union_of_hierarchy: Dict[str, Set[str]],
    true_taxon_lineage: Set[str],
    lengths: List[int],
    num_samples: int,
    threshold: float,
    seed: int,
) -> List[dict]:
    """Run the full sampling experiment across all contig lengths."""
    rng = random.Random(seed)
    prediction_entry_names = set(predictions.keys())
    results = []

    for length in lengths:
        if length > genome_size:
            print(f"[WARN] Contig length {length} exceeds genome size {genome_size}, skipping")
            continue

        print(f"\n{'='*60}")
        print(f"Contig length: {length:,} bp  ({num_samples} samples)")
        print(f"{'='*60}")

        for sample_idx in range(num_samples):
            contig_start, protein_names = sample_contig_proteins(
                cds_entries, genome_size, length, rng
            )

            proteins_with_preds = [p for p in protein_names if p in prediction_entry_names]

            if len(proteins_with_preds) == 0:
                results.append({
                    'length': length,
                    'sample_idx': sample_idx,
                    'contig_start': contig_start,
                    'num_proteins': 0,
                    'num_proteins_with_preds': 0,
                    'is_matching': None,
                    'most_specific_taxon': None,
                    'num_predicted_taxons': 0,
                    'total_flips': 0,
                    'status': 'no_proteins',
                })
                continue

            contig_predictions = {p: predictions[p] for p in proteins_with_preds}

            annotated_predictions = get_annotated_predictions(contig_predictions, threshold)
            num_annotated = sum(len(terms) for terms in annotated_predictions.values())
            if num_annotated == 0:
                results.append({
                    'length': length,
                    'sample_idx': sample_idx,
                    'contig_start': contig_start,
                    'num_proteins': len(protein_names),
                    'num_proteins_with_preds': len(proteins_with_preds),
                    'is_matching': None,
                    'most_specific_taxon': None,
                    'num_predicted_taxons': 0,
                    'total_flips': 0,
                    'status': 'no_annotations_above_threshold',
                })
                continue

            flip_cost_per_go_term, _ = compute_demotion_flip_cost(contig_predictions, threshold)

            try:
                adjusted_preds, predicted_taxons, total_flips = solve_sat_ortools_taxon_assignment(
                    annotated_predictions,
                    in_taxon_constraints,
                    never_in_taxon_constraints,
                    go_hierarchy,
                    flip_cost_per_go_term,
                    taxon_is_a_hierarchy,
                    taxon_disjoint_from_hierarchy,
                    taxon_union_of_hierarchy,
                    genome_taxon_id=None,
                )
            except Exception as e:
                print(f"  [ERROR] Sample {sample_idx}: solver failed: {e}")
                results.append({
                    'length': length,
                    'sample_idx': sample_idx,
                    'contig_start': contig_start,
                    'num_proteins': len(protein_names),
                    'num_proteins_with_preds': len(proteins_with_preds),
                    'is_matching': None,
                    'most_specific_taxon': None,
                    'num_predicted_taxons': 0,
                    'total_flips': 0,
                    'status': 'solver_error',
                })
                continue

            matching = is_matching_taxon(predicted_taxons, true_taxon_lineage)
            most_specific = get_most_specific_taxon(predicted_taxons, taxon_is_a_hierarchy)

            results.append({
                'length': length,
                'sample_idx': sample_idx,
                'contig_start': contig_start,
                'num_proteins': len(protein_names),
                'num_proteins_with_preds': len(proteins_with_preds),
                'is_matching': matching,
                'most_specific_taxon': most_specific,
                'num_predicted_taxons': len(predicted_taxons),
                'total_flips': total_flips,
                'status': 'ok',
            })

            status_str = "MATCH" if matching else "MISMATCH"
            print(f"  Sample {sample_idx:3d}: {len(proteins_with_preds):4d} proteins, "
                  f"{status_str}, most_specific={most_specific}, "
                  f"flips={total_flips}, #taxons={len(predicted_taxons)}")

    return results


def load_detail_tsv(detail_path: str) -> List[dict]:
    """Load a previously saved detail TSV back into a list of result dicts."""
    results = []
    with open(detail_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            is_matching_raw = row.get('is_matching', '')
            if is_matching_raw == 'True':
                is_matching = True
            elif is_matching_raw == 'False':
                is_matching = False
            else:
                is_matching = None
            most_specific = row.get('most_specific_taxon', '') or None
            results.append({
                'length': int(row['length']),
                'sample_idx': int(row['sample_idx']),
                'contig_start': int(row['contig_start']),
                'num_proteins': int(row['num_proteins']),
                'num_proteins_with_preds': int(row['num_proteins_with_preds']),
                'is_matching': is_matching,
                'most_specific_taxon': most_specific,
                'num_predicted_taxons': int(row['num_predicted_taxons']),
                'total_flips': int(row['total_flips']),
                'status': row['status'],
            })
    print(f"[INFO] Loaded {len(results)} rows from {detail_path}")
    return results


def load_genome_size_from_summary(summary_path: str) -> Optional[int]:
    """Load genome_size from an existing summary TSV if present."""
    if not os.path.isfile(summary_path):
        return None
    with open(summary_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        row = next(reader, None)
        if row is None or 'genome_size' not in row:
            return None
        try:
            return int(row['genome_size'])
        except (ValueError, KeyError):
            return None


def save_results(
    results: List[dict],
    output_dir: str,
    taxon_id: str,
    suffix: str = '',
    genome_size: Optional[int] = None,
):
    """Save per-sample results and summary TSV."""
    os.makedirs(output_dir, exist_ok=True)

    detail_path = os.path.join(output_dir, f'metagenomics_detail_{taxon_id}{suffix}.tsv')
    with open(detail_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'length', 'sample_idx', 'contig_start', 'num_proteins',
            'num_proteins_with_preds', 'is_matching', 'most_specific_taxon',
            'num_predicted_taxons', 'total_flips', 'status',
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    print(f"[INFO] Wrote detailed results to {detail_path}")

    summary = compute_summary(results, genome_size=genome_size)
    summary_fieldnames = [
        'length', 'total_samples', 'valid_samples', 'matching',
        'accuracy', 'general_matching', 'general_accuracy',
        'avg_proteins', 'avg_proteins_with_preds', 'avg_flips',
    ]
    if genome_size is not None:
        summary_fieldnames.append('genome_size')
    summary_path = os.path.join(output_dir, f'metagenomics_summary_{taxon_id}{suffix}.tsv')
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=summary_fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary)
    print(f"[INFO] Wrote summary results to {summary_path}")

    return summary


def compute_summary(results: List[dict], genome_size: Optional[int] = None) -> List[dict]:
    """Aggregate per-sample results into per-length summary."""
    by_length = defaultdict(list)
    for r in results:
        by_length[r['length']].append(r)

    summary = []
    for length in sorted(by_length.keys()):
        samples = by_length[length]
        valid = [s for s in samples if s['status'] == 'ok']
        matching = sum(1 for s in valid if s['is_matching'])
        accuracy = matching / len(valid) if valid else 0.0
        general_matching = sum(
            1 for s in valid
            if is_general_matching(s['is_matching'], s.get('most_specific_taxon'))
        )
        general_accuracy = general_matching / len(valid) if valid else 0.0
        avg_proteins = (sum(s['num_proteins'] for s in samples) / len(samples)) if samples else 0
        avg_proteins_with_preds = (
            sum(s['num_proteins_with_preds'] for s in samples) / len(samples)
        ) if samples else 0
        avg_flips = (sum(s['total_flips'] for s in samples) / len(samples)) if samples else 0
        row = {
            'length': length,
            'total_samples': len(samples),
            'valid_samples': len(valid),
            'matching': matching,
            'accuracy': f"{accuracy:.4f}",
            'general_matching': general_matching,
            'general_accuracy': f"{general_accuracy:.4f}",
            'avg_proteins': f"{avg_proteins:.1f}",
            'avg_proteins_with_preds': f"{avg_proteins_with_preds:.1f}",
            'avg_flips': f"{avg_flips:.0f}",
        }
        if genome_size is not None:
            row['genome_size'] = genome_size
        summary.append(row)
    return summary


def _draw_summary_axes(ax1, summary: List[dict]):
    """Draw accuracy bars and avg-proteins line onto a pair of axes (ax1, ax1.twinx())."""
    lengths = [int(s['length']) for s in summary]
    accuracies = [float(s['accuracy']) for s in summary]
    general_accuracies = [float(s.get('general_accuracy', 0)) for s in summary]
    valid_counts = [int(s['valid_samples']) for s in summary]
    avg_prots = [float(s['avg_proteins_with_preds']) for s in summary]

    # General matching is always a subset of strict matching in practice.
    # Blue (bottom): specific correct matches (strict but not general).
    # Green (top, stacked): trivially/vaguely correct matches (general = NCBITaxon_131567 or empty).
    # Total bar height = strict accuracy.
    green_heights = general_accuracies
    blue_heights = [max(0.0, a - g) for a, g in zip(accuracies, general_accuracies)]
    bar_tops = accuracies

    x_labels = [f"{l:,}" for l in lengths]
    x_pos = list(range(len(lengths)))

    bar_width = 0.6
    ax1.bar(x_pos, blue_heights, width=bar_width, color='steelblue', alpha=0.85,
            label='Specific matches (is_matching, specific taxon)')
    ax1.bar(x_pos, green_heights, width=bar_width, bottom=blue_heights,
            color='mediumseagreen', alpha=0.85,
            label='General matches (NCBITaxon_131567 or no specific taxon)')

    ax1.set_xlabel('Contig Length (bp)')
    ax1.set_ylabel('Taxon Matching Accuracy')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(x_labels, rotation=15)
    ax1.set_ylim(0, 1.15)
    ax1.legend(loc='upper left', fontsize=8)

    for x, top, acc, gen_acc, n_valid in zip(x_pos, bar_tops, accuracies, general_accuracies, valid_counts):
        ax1.text(x, top + 0.02,
                 f"acc={acc:.2f}\ngen={gen_acc:.2f}\n(n={n_valid})", ha='center', va='bottom', fontsize=8,
                 color='steelblue')

    ax2 = ax1.twinx()
    ax2.plot(x_pos, avg_prots, 'o-', color='darkorange', label='Avg proteins')
    ax2.set_ylabel('Avg Proteins with Predictions', color='darkorange')
    ax2.tick_params(axis='y', labelcolor='darkorange')
    for x, avg in zip(x_pos, avg_prots):
        ax2.annotate(
            f"{avg:.1f}",
            (x, avg),
            textcoords="offset points",
            xytext=(0, -12),
            ha='center',
            va='top',
            color='darkorange',
            fontsize=9,
        )


def plot_results(
    summary: List[dict],
    output_dir: str,
    taxon_id: str,
    suffix: str = '',
    genome_size: Optional[int] = None,
):
    """Create bar chart of accuracy vs contig length."""
    fig, ax1 = plt.subplots(figsize=(10, 6))
    _draw_summary_axes(ax1, summary)
    title = f'Metagenomics Taxon Prediction Accuracy\nTaxon {taxon_id}'
    if genome_size is not None:
        title += f' (genome: {genome_size:,} bp)'
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()

    plot_path = os.path.join(output_dir, f'metagenomics_accuracy_{taxon_id}{suffix}.png')
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"[INFO] Saved plot to {plot_path}")


def plot_results_combined(
    summaries_by_onto: Dict[str, List[dict]],
    output_dir: str,
    taxon_id: str,
    suffix: str = '',
    genome_size: Optional[int] = None,
):
    """Create a single figure with one subplot per subontology."""
    ontos = [o for o in SUBONTOLOGIES if o in summaries_by_onto]
    if not ontos:
        return

    fig, axes = plt.subplots(1, len(ontos), figsize=(10 * len(ontos), 6), squeeze=False)
    for col, onto in enumerate(ontos):
        ax1 = axes[0][col]
        _draw_summary_axes(ax1, summaries_by_onto[onto])
        ax1.set_title(onto.upper(), fontsize=12, fontweight='bold')

    title = f'Metagenomics Taxon Prediction Accuracy — Taxon {taxon_id}'
    if genome_size is not None:
        title += f' (genome: {genome_size:,} bp)'
    fig.suptitle(title, fontsize=14)
    fig.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    plot_path = os.path.join(output_dir, f'metagenomics_accuracy_{taxon_id}{suffix}_combined.png')
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"[INFO] Saved combined plot to {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Sample contigs and evaluate taxon prediction accuracy.'
    )
    parser.add_argument('--predictions-file', default=None,
                        help='Path to predictions TSV (protein_entry_name<TAB>GO:term|score...)')
    parser.add_argument('--predictions-dir', default=None,
                        help='Directory containing subontology subdirectories with predictions TSV files. When given, predictions-file is ignored.')
    parser.add_argument('--cds-file', required=True,
                        help='Path to CDS positions TSV (start, end, complement, entry_name, is_swissprot)')
    parser.add_argument('--taxon-id', required=True,
                        help='True taxon ID (e.g., 83332)')
    parser.add_argument('--constraints', default=DEFAULT_CONSTRAINTS,
                        help='Path to GO taxon constraints file')
    parser.add_argument('--go-hierarchy', default=DEFAULT_GO_HIERARCHY,
                        help='Path to GO hierarchy file')
    parser.add_argument('--taxon-hierarchy', default=DEFAULT_TAXON_HIERARCHY,
                        help='Path to taxon hierarchy file')
    parser.add_argument('--ncbitaxon-hierarchy', default=DEFAULT_NCBITAXON_HIERARCHY,
                        help='Path to NCBITaxon hierarchy file')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Prediction score threshold (default: 0.5)')
    parser.add_argument('--num-samples', type=int, default=NUM_SAMPLES,
                        help=f'Number of samples per contig length (default: {NUM_SAMPLES})')
    parser.add_argument('--lengths', type=str, default=','.join(str(l) for l in LENGTHS),
                        help=f'Comma-separated contig lengths (default: {",".join(str(l) for l in LENGTHS)})')
    parser.add_argument('--output-dir', default='metagenomics_results',
                        help='Output directory (default: metagenomics_results)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('--replot', default=None,
                        help='Path to an existing detail TSV (or directory containing detail TSVs '
                             'named metagenomics_detail_<taxon_id>.tsv per subontology). '
                             'Skips the experiment and re-generates plots from saved data.')
    parser.add_argument('--full-sequence-only', action='store_true', default=False,
                        help='Run the solver once on the entire genome sequence (all proteins) '
                             'instead of sampling contigs. Overrides --lengths and --num-samples.')
    args = parser.parse_args()

    lengths = [int(x.strip()) for x in args.lengths.split(',')]

    # --- Replot-only mode ---
    if args.replot is not None:
        replot_path = args.replot
        if os.path.isfile(replot_path):
            results = load_detail_tsv(replot_path)
            summary_dir = os.path.dirname(replot_path)
            detail_basename = os.path.basename(replot_path)
            # Derive summary path: metagenomics_detail_X.tsv -> metagenomics_summary_X.tsv
            summary_basename = detail_basename.replace('metagenomics_detail_', 'metagenomics_summary_')
            summary_file = os.path.join(summary_dir, summary_basename)
            genome_size = load_genome_size_from_summary(summary_file)
            summary = compute_summary(results, genome_size=genome_size)
            os.makedirs(args.output_dir, exist_ok=True)
            save_results(results, args.output_dir, args.taxon_id, genome_size=genome_size)
            plot_results(summary, args.output_dir, args.taxon_id, genome_size=genome_size)
        elif os.path.isdir(replot_path):
            summaries_by_onto: Dict[str, List[dict]] = {}
            genome_size: Optional[int] = None
            for onto in SUBONTOLOGIES:
                detail_file = os.path.join(
                    replot_path, onto,
                    f'metagenomics_detail_{args.taxon_id}.tsv'
                )
                if not os.path.isfile(detail_file):
                    print(f"[WARN] No detail TSV found for '{onto}' at {detail_file}, skipping")
                    continue
                if genome_size is None:
                    # Try summary with or without _full suffix
                    for suf in ['', '_full']:
                        summary_file = os.path.join(
                            replot_path, onto,
                            f'metagenomics_summary_{args.taxon_id}{suf}.tsv'
                        )
                        genome_size = load_genome_size_from_summary(summary_file)
                        if genome_size is not None:
                            break
                results = load_detail_tsv(detail_file)
                onto_output_dir = os.path.join(args.output_dir, onto)
                summaries_by_onto[onto] = compute_summary(results, genome_size=genome_size)
                save_results(results, onto_output_dir, args.taxon_id, genome_size=genome_size)
                plot_results(
                    summaries_by_onto[onto], onto_output_dir, args.taxon_id,
                    genome_size=genome_size,
                )
            plot_results_combined(
                summaries_by_onto, args.output_dir, args.taxon_id,
                genome_size=genome_size,
            )
        else:
            print(f"[ERROR] --replot path does not exist: {replot_path}")
            sys.exit(1)
        print("\nDone!")
        return

    if args.predictions_dir is None and args.predictions_file is None:
        print("[ERROR] Either --predictions-file or --predictions-dir must be given.")
        sys.exit(1)

    # --- Load shared data ---
    print("=" * 60)
    print("Loading data...")
    print("=" * 60)

    print(f"\n[1/5] Loading CDS positions from {args.cds_file}...")
    cds_entries = load_cds_positions(args.cds_file)
    if not cds_entries:
        print("[ERROR] No CDS entries loaded. Check that the CDS file is valid and non-empty.")
        sys.exit(1)
    genome_size = max(end for _, end, _ in cds_entries)
    print(f"  Genome size (from CDS): {genome_size:,} bp")

    if args.full_sequence_only:
        print(f"\n[INFO] --full-sequence-only: running solver once on all {len(cds_entries)} CDS entries "
              f"(genome_size={genome_size:,} bp), ignoring --lengths and --num-samples.")
        lengths = [genome_size]
        args.num_samples = 1

    print(f"\n[3/5] Loading GO taxon constraints from {args.constraints}...")
    in_taxon_constraints, never_in_taxon_constraints = load_constraints(args.constraints)
    print(f"  {len(in_taxon_constraints)} only_in_taxon, {len(never_in_taxon_constraints)} never_in_taxon")

    print(f"\n[4/5] Loading GO hierarchy from {args.go_hierarchy}...")
    go_hierarchy = load_go_hierarchy(args.go_hierarchy) if args.go_hierarchy else defaultdict(set)
    print(f"  {len(go_hierarchy)} GO hierarchy relationships")

    print(f"\n[5/5] Loading taxon hierarchy from {args.taxon_hierarchy}...")
    taxon_is_a_hierarchy, taxon_disjoint_from_hierarchy, taxon_union_of_hierarchy = (
        load_taxon_hierarchy(args.taxon_hierarchy, args.ncbitaxon_hierarchy)
        if args.taxon_hierarchy else (defaultdict(set), defaultdict(set), defaultdict(set))
    )
    print(f"  {len(taxon_is_a_hierarchy)} is_a, {len(taxon_disjoint_from_hierarchy)} disjoint_from, "
          f"{len(taxon_union_of_hierarchy)} union_of")

    normalized_taxon = normalize_taxon_id(args.taxon_id)
    true_taxon_lineage = get_all_ancestors(normalized_taxon, taxon_is_a_hierarchy)
    print(f"\n  True taxon lineage for {args.taxon_id} ({normalized_taxon}): {len(true_taxon_lineage)} ancestors")

    def _run_for_predictions(pred_file: str, output_dir: str, save_individual_plot: bool = True) -> List[dict]:
        print(f"\n[2/5] Loading predictions from {pred_file}...")
        predictions = load_predictions(pred_file)
        print(f"  Loaded predictions for {len(predictions)} proteins")
        matched_proteins = set(e[2] for e in cds_entries) & set(predictions.keys())
        print(f"  CDS proteins with predictions: {len(matched_proteins)}/{len(cds_entries)}")

        print("\n" + "=" * 60)
        print("Running metagenomics contig experiment...")
        print(f"  Lengths: {lengths}")
        print(f"  Samples per length: {args.num_samples}")
        print(f"  Threshold: {args.threshold}")
        print(f"  Seed: {args.seed}")
        print("=" * 60)

        results = run_metagenomics_experiment(
            predictions=predictions,
            cds_entries=cds_entries,
            genome_size=genome_size,
            in_taxon_constraints=in_taxon_constraints,
            never_in_taxon_constraints=never_in_taxon_constraints,
            go_hierarchy=go_hierarchy,
            taxon_is_a_hierarchy=taxon_is_a_hierarchy,
            taxon_disjoint_from_hierarchy=taxon_disjoint_from_hierarchy,
            taxon_union_of_hierarchy=taxon_union_of_hierarchy,
            true_taxon_lineage=true_taxon_lineage,
            lengths=lengths,
            num_samples=args.num_samples,
            threshold=args.threshold,
            seed=args.seed,
        )

        file_suffix = '_full' if args.full_sequence_only else ''
        summary = save_results(
            results, output_dir, args.taxon_id, suffix=file_suffix, genome_size=genome_size
        )
        if save_individual_plot:
            plot_results(
                summary, output_dir, args.taxon_id, suffix=file_suffix, genome_size=genome_size
            )

        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        for s in summary:
            print(f"  Length {int(s['length']):>10,}: accuracy={s['accuracy']}, "
                  f"valid={s['valid_samples']}/{s['total_samples']}, "
                  f"avg_proteins={s['avg_proteins_with_preds']}, "
                  f"avg_flips={s['avg_flips']}")

        return summary

    # --- Run experiment(s) ---
    if args.predictions_dir:
        summaries_by_onto: Dict[str, List[dict]] = {}
        for onto in SUBONTOLOGIES:
            pred_file = find_predictions_file(args.predictions_dir, onto, args.taxon_id)
            if pred_file is None:
                print(f"\n[WARN] No predictions file found for subontology '{onto}' "
                      f"under {args.predictions_dir}/{onto}/, skipping")
                continue
            print(f"\n{'#'*60}")
            print(f"# Subontology: {onto.upper()}")
            print(f"{'#'*60}")
            onto_output_dir = os.path.join(args.output_dir, onto)
            summaries_by_onto[onto] = _run_for_predictions(
                pred_file, onto_output_dir, save_individual_plot=False
            )
        plot_results_combined(
            summaries_by_onto, args.output_dir, args.taxon_id,
            suffix='_full' if args.full_sequence_only else '',
            genome_size=genome_size,
        )
    else:
        _run_for_predictions(args.predictions_file, args.output_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
