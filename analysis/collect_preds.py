"""
This script collects the predictions for the test proteins from given prediction files.
Computes CAFA-like metrics.

The predictions are collected in the format: protein_id<tab>GO:term|score<tab>...
"""

from typing import List, Tuple, Dict, Optional
from collections import defaultdict
import os
import glob
import json
import sys
import argparse

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PIPELINE_DIR = os.path.join(BASE_DIR, 'pipeline')
GAEF_DIR = os.path.abspath(os.path.join(BASE_DIR, '..', 'GAEF'))

def load_test_proteins(test_proteins_file: str) -> Dict[str, List[str]]:
    """
    Load the test proteins from the test proteins file with matching taxon IDs.

    Example: (proteins_by_date_23-MAY-2024.tsv file)
    protein_id	protein_name	integration_date	taxon_id	taxon_id_matched
    10H_STRNX	Strychnine-10-hydroxylase {ECO:0000303|PubMed:35794473}	27-NOV-2024	NA	False
    34K1_AEDAL	Salivary protein SG34 {ECO:0000305}	02-OCT-2024	7160	True
    ...

    Example:
    protein_id taxon_id
    10H_STRNX NA
    34K1_AEDAL 7160
    ...

    Returns:
    Dict[str, List[str]]: Dictionary mapping taxon IDs to list of test protein IDs.
    """
    test_proteins = defaultdict(list)
    with open(test_proteins_file, 'r') as f:
        # obtain taxon_id index from the header
        header = f.readline().strip().split('\t')
        taxon_id_index = header.index('taxon_id')
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            protein_id = parts[0]
            if protein_id == 'protein_id': # skip header
                continue
            taxon_id = parts[taxon_id_index]
            if taxon_id != 'NA':
                test_proteins[taxon_id].append(protein_id)
    return test_proteins

def load_predictions(predictions_dir: str, test_proteins: Dict[str, List[str]]) -> Dict[str, Dict[str, float]]:
    """
    Load the predictions from the predictions directory for the test proteins.
    The predictions are in the format: protein_id<tab>GO:term|score<tab>...
    Returns:
    Dict[str, Dict[str, float]]: Dictionary mapping protein_id to its predictions.
    """
    predictions = defaultdict(dict)

    for taxon_id, protein_ids in test_proteins.items():
        # find predictions file for the given taxon_id
        predictions_files = glob.glob(os.path.join(predictions_dir, f'*_taxon_{taxon_id}.tsv'))
        if len(predictions_files) == 0:
            print(f"Warning: Predictions file for taxon {taxon_id} not found")
            continue
        predictions_file = predictions_files[0]

        protein_id_set = set(protein_ids)
        with open(predictions_file, 'r') as f:
            for line in f:
                protein_id = line.split('\t', 1)[0]
                if protein_id in protein_id_set:
                    _, *preds = line.strip().split('\t')  # obtain the predictions for the protein
                    predictions[protein_id] = {
                        go_id: float(score_str)
                        for pred in preds
                        for go_id, score_str in [pred.split('|', 1)]
                    }

    return predictions


def load_predictions_from_single_file(
    predictions_file: str,
    test_proteins: Dict[str, List[str]],
) -> Dict[str, Dict[str, float]]:
    """
    Load predictions from a single ProtGO output file (e.g. protgo_predictions_filtered_mf.tsv).
    Format: protein_id<tab>GO:term|score<tab>...
    Returns:
    Dict[str, Dict[str, float]]: Dictionary mapping protein_id to its predictions.
    """
    protein_id_set = set()
    for protein_ids in test_proteins.values():
        protein_id_set.update(protein_ids)

    predictions = {}
    if not os.path.exists(predictions_file):
        print(f"Warning: Predictions file not found: {predictions_file}")
        return predictions

    with open(predictions_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            protein_id = parts[0]
            if protein_id in protein_id_set:
                preds = parts[1:]
                predictions[protein_id] = {
                    go_id: float(score_str)
                    for pred in preds
                    for go_id, score_str in [pred.split('|', 1)]
                }

    return predictions


def compute_cafa_metrics(
    predictions: Dict[str, Dict[str, float]],
    test_proteins: Dict[str, List[str]],
    annotations_dir: str,
    go_file: str,
    subontology: str,
) -> Dict:
    """
    Compute CAFA-like metrics using pooled predictions/annotations.
    """
    if PIPELINE_DIR not in sys.path:
        sys.path.append(PIPELINE_DIR)
    if BASE_DIR not in sys.path:
        sys.path.append(BASE_DIR)
    if GAEF_DIR not in sys.path:
        sys.path.append(GAEF_DIR)

    import importlib.util

    genome_scores_path = os.path.join(PIPELINE_DIR, 'genome_scores_evaluator.py')
    spec = importlib.util.spec_from_file_location("genome_scores_evaluator", genome_scores_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load genome_scores_evaluator from {genome_scores_path}")
    genome_scores = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(genome_scores)

    load_annotations = genome_scores.load_annotations
    evaluate_prediction_metrics = genome_scores.evaluate_prediction_metrics
    Ontology = genome_scores.Ontology

    if not os.path.exists(annotations_dir):
        print(f"[Error] Annotations directory not found: {annotations_dir}")
        return {}
    if not os.path.exists(go_file):
        print(f"[Error] GO file not found: {go_file}")
        return {}

    # Build protein -> taxon mapping to avoid ID collisions
    protein_to_taxon = {}
    for taxon_id, protein_ids in test_proteins.items():
        for protein_id in protein_ids:
            protein_to_taxon[protein_id] = taxon_id

    # Pool predictions with taxon prefix
    pooled_predictions: Dict[str, Dict[str, float]] = {}
    for protein_id, scores in predictions.items():
        taxon_id = protein_to_taxon.get(protein_id)
        if taxon_id is None:
            continue
        pooled_id = f"{taxon_id}:{protein_id}"
        pooled_predictions[pooled_id] = scores

    # Pool annotations with taxon prefix
    pooled_annotations: Dict[str, set] = {}
    missing_taxa = []
    for taxon_id, protein_ids in test_proteins.items():
        annotation_file = os.path.join(annotations_dir, f'annots_taxon_{taxon_id}.tsv')
        if not os.path.exists(annotation_file):
            missing_taxa.append(taxon_id)
            continue
        taxon_annotations = load_annotations(annotation_file)
        protein_set = set(protein_ids)
        for prot_id, go_terms in taxon_annotations.items():
            if prot_id in protein_set:
                pooled_id = f"{taxon_id}:{prot_id}"
                pooled_annotations[pooled_id] = go_terms

    if not pooled_annotations:
        print("[Error] No annotations found for provided test proteins")
        return {}

    go = Ontology(go_file)

    metrics_by_ont = {}
    metrics = evaluate_prediction_metrics(go, subontology, pooled_annotations, pooled_predictions)
    if metrics is not None:
        metrics_by_ont[subontology] = metrics

    summary = {}
    metric_keys = ['fmax', 'smin', 'avg_auc', 'aupr', 'precision', 'recall', 'sensitivity', 'specificity', 'wfmax', 'avgic']
    for key in metric_keys:
        values = [m[key] for m in metrics_by_ont.values() if key in m]
        if values:
            summary[f'mean_{key}'] = sum(values) / len(values)

    return {
        'per_ontology': metrics_by_ont,
        'summary': summary,
        'total_annotated_proteins': len(pooled_annotations),
        'missing_annotation_taxa': missing_taxa,
    }


def main(
    predictions_dir: Optional[str],
    predictions_file: Optional[str],
    test_proteins_file: str,
    output_file: str,
    annotations_dir: str,
    go_file: str,
    subontology: str,
):
    """
    Main function to collect predictions and compute CAFA-like metrics.
    """
    if predictions_file and predictions_dir:
        raise ValueError("Provide either --predictions_file or --predictions_dir, not both")
    if not predictions_file and not predictions_dir:
        raise ValueError("Provide either --predictions_file or --predictions_dir")

    test_proteins = load_test_proteins(test_proteins_file)

    if predictions_file:
        predictions = load_predictions_from_single_file(predictions_file, test_proteins)
        print("=== Predictions file: ", predictions_file)
    else:
        predictions = load_predictions(predictions_dir, test_proteins)
        print("=== Predictions directory: ", predictions_dir)
    print(f"Loaded {len(predictions)} proteins for {len(test_proteins)} test organisms (including not found for some taxa)")

    cafa_metrics = compute_cafa_metrics(predictions, test_proteins, annotations_dir, go_file, subontology)
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(cafa_metrics, f, indent=2)

    print("=== Micro-Average Metrics ===")
    for ontology, metrics in cafa_metrics['per_ontology'].items():
        print(f"{ontology}:")
        print(f"  Fmax: {metrics['fmax']:.5f}")
        print(f"  Smin: {metrics['smin']:.5f}")
        print(f"  AUC: {metrics['avg_auc']:.5f}")
        print(f"  AUPR: {metrics['aupr']:.5f}")
        print(f"  Precision: {metrics['precision']:.5f}")
        print(f"  Recall: {metrics['recall']:.5f}")
        print(f"  Sensitivity (term-centric): {metrics.get('sensitivity', 0.0):.5f}")
        print(f"  Specificity (term-centric): {metrics.get('specificity', 0.0):.5f}")
        print(f"  TP: {metrics.get('tp', 'N/A')}, FP: {metrics.get('fp', 'N/A')}, FN: {metrics.get('fn', 'N/A')}")
        print(f"  Threshold: {metrics['tmax']}")
        print(f"  WFmax: {metrics['wfmax']:.5f}")
        print(f"  AVGIC: {metrics['avgic']:.5f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Collect predictions and compute CAFA-like metrics')
    parser.add_argument('--predictions_dir', type=str, default=None, help='Directory containing per-taxon predictions (*_taxon_{taxon_id}.tsv)')
    parser.add_argument('--predictions_file', type=str, default=None, help='Single ProtGO predictions file (e.g. protgo_predictions_filtered_mf.tsv)')
    parser.add_argument('--test_proteins_file', type=str, default='proteins_by_date_23-MAY-2024.tsv', help='File containing test proteins')
    parser.add_argument('--output_file', type=str, default='micro_avg_metrics.json', help='File to save the metrics')
    parser.add_argument('--annotations_dir', type=str, default="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic", help='Directory containing annotations')
    parser.add_argument('--go_file', type=str, default="data/go-basic.obo", help='Path to GO file')
    parser.add_argument('--subontology', type=str, default='cc', choices=['cc', 'bp', 'mf'], help='Subontology to evaluate')
    args = parser.parse_args()
    main(args.predictions_dir, args.predictions_file, args.test_proteins_file, args.output_file, args.annotations_dir, args.go_file, args.subontology)