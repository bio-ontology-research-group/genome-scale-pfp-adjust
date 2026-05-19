"""
Given the test predictions and optimized annotations, obtain the adjusted proteins (i.e. proteins that have different predictions in the test and optimized annotations).

Output a tsv file with protein_id and taxon_id columns.

"""

import os
import argparse
import csv
from typing import Dict
from collections import defaultdict
from tqdm import tqdm


def parse_prediction_file(filepath: str) -> Dict[str, Dict[str, str]]:
    """
    Parse a TSV prediction file into a mapping of protein_id -> dict of GO term -> score.

    Each line has the format:
        protein_id\tGO:XXXXXXX|score\tGO:YYYYYYY|score\t...
    """
    result: Dict[str, Dict[str, str]] = {}
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row:
                continue
            protein_id = row[0]
            go_terms: Dict[str, str] = {}
            for field in row[1:]:
                if '|' in field:
                    go_term, score = field.split('|', 1)
                    go_terms[go_term] = score
            result[protein_id] = go_terms
    return result


def extract_taxon_id(filename: str) -> str:
    """
    Extract taxon_id from a prediction/optimized filename.

    Expected formats:
        predictions_fold_01_taxon_83332.tsv
        optimized_fold_01_taxon_83332.tsv
    """
    # Strip extension
    base = os.path.splitext(filename)[0]
    # taxon_id follows '_taxon_'
    if '_taxon_' not in base:
        raise ValueError(f"Cannot extract taxon_id from filename: {filename}")
    return base.split('_taxon_')[1]


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Compare original test predictions with adjusted optimized predictions '
            'and output proteins whose GO term sets have changed.'
        )
    )
    parser.add_argument(
        '--predictions-dir', default="${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results/cc/predictions/",
        help='Directory containing original prediction TSV files '
             '(predictions_fold_XX_taxon_YYYYYYY.tsv).'
    )
    parser.add_argument(
        '--optimized-dir', default=None,
        help='Directory containing adjusted optimized TSV files '
             '(optimized_fold_XX_taxon_YYYYYYY.tsv).'
    )
    parser.add_argument(
        '--annotations-dir', default="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic",
        help='Directory containing per-taxon annotation files.'
    )
    parser.add_argument(
        '--output-file', default=None,
        help='Output TSV file path (columns: protein_id, taxon_id).'
    )
    args = parser.parse_args()

    if not os.path.isdir(args.predictions_dir):
        raise FileNotFoundError(f"predictions-dir not found: {args.predictions_dir}")
    if args.optimized_dir is None:
        args.optimized_dir = args.predictions_dir.replace('predictions/', 'optimized/')
    if not os.path.isdir(args.optimized_dir):
        raise FileNotFoundError(f"optimized-dir not found: {args.optimized_dir}")

    if args.output_file is None:
        args.output_file = os.path.join(args.predictions_dir, '..', 'adjusted_proteins.tsv')
    os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)

    adjusted_proteins = []  # list of (protein_id, taxon_id)
    total_proteins = 0

    optimized_files = [
        f for f in os.listdir(args.optimized_dir)
        if f.startswith('optimized_') and f.endswith('.tsv')
    ]

    print(f"[INFO] Found {len(optimized_files)} optimized files in {args.optimized_dir}")

    for opt_filename in tqdm(sorted(optimized_files), desc="Processing optimized files"):
        taxon_id = extract_taxon_id(opt_filename)

        # Derive the matching predictions filename
        pred_filename = opt_filename.replace('optimized_', 'predictions_', 1)
        opt_path = os.path.join(args.optimized_dir, opt_filename)
        pred_path = os.path.join(args.predictions_dir, pred_filename)

        if not os.path.exists(pred_path):
            print(f"[WARNING] No matching predictions file for {opt_filename}: expected {pred_path}")
            continue

        optimized = parse_prediction_file(opt_path)
        predictions = parse_prediction_file(pred_path)

        total_proteins += len(predictions)

        for protein_id in predictions:
            opt_dict = optimized.get(protein_id, {})
            pred_dict = predictions.get(protein_id, {})
            opt_terms = set(opt_dict.keys())
            pred_terms = set(pred_dict.keys())
            if opt_terms != pred_terms:
                removed_terms = pred_terms - opt_terms
                removed_str = ",".join([f"{t}|{pred_dict[t]}" for t in removed_terms])
                adjusted_proteins.append((protein_id, taxon_id, removed_str))

    print(f"[INFO] Total adjusted proteins: {len(adjusted_proteins)} out of {total_proteins} total proteins ({(len(adjusted_proteins) / total_proteins) * 100:.2f}%)")

    taxon_to_proteins = defaultdict(list)
    for protein_id, taxon_id, removed_str in adjusted_proteins:
        taxon_to_proteins[taxon_id].append((protein_id, removed_str))

    filtered_adjusted_proteins = []

    for taxon_id, protein_ids in taxon_to_proteins.items():
        annot_file = os.path.join(args.annotations_dir, f"annots_taxon_{taxon_id}.tsv")
        valid_proteins = set()
        if os.path.exists(annot_file):
            with open(annot_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) > 1 and any(t for t in parts[1:] if t):
                        valid_proteins.add(parts[0])

        for protein_id, removed_str in protein_ids:
            if protein_id in valid_proteins:
                filtered_adjusted_proteins.append((protein_id, taxon_id, removed_str))

    print(f"[INFO] Total adjusted proteins after filtering: {len(filtered_adjusted_proteins)} (removed {len(adjusted_proteins) - len(filtered_adjusted_proteins)})")

    with open(args.output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(['protein_id', 'taxon_id', 'removed_terms'])
        for protein_id, taxon_id, removed_str in filtered_adjusted_proteins:
            writer.writerow([protein_id, taxon_id, removed_str])

    print(f"[INFO] Written to {args.output_file}")


if __name__ == '__main__':
    main()
