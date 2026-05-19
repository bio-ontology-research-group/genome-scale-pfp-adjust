"""
Compare term-specific sensitivity and specificity between predictions and optimized annotations.
Reports top terms by largest absolute difference, with ground-truth annotation counts.
"""

from typing import Dict, List, Optional
from collections import defaultdict
import os
import glob
import sys
import argparse
import math
import statistics

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
PIPELINE_DIR = os.path.join(BASE_DIR, 'pipeline')
GAEF_DIR = os.path.abspath(os.path.join(BASE_DIR, '..', 'GAEF'))


def load_test_proteins(test_proteins_file: str) -> Dict[str, List[str]]:
    """Load test proteins from TSV with taxon_id column."""
    test_proteins = defaultdict(list)
    with open(test_proteins_file, 'r') as f:
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
            if protein_id == 'protein_id':
                continue
            taxon_id = parts[taxon_id_index]
            if taxon_id != 'NA':
                test_proteins[taxon_id].append(protein_id)
    return test_proteins


def load_predictions_from_dir(
    predictions_dir: str, test_proteins: Dict[str, List[str]]
) -> Dict[str, Dict[str, float]]:
    """Load predictions from per-taxon TSV files."""
    predictions = defaultdict(dict)
    for taxon_id, protein_ids in test_proteins.items():
        files = glob.glob(os.path.join(predictions_dir, f'*_taxon_{taxon_id}.tsv'))
        if not files:
            continue
        protein_id_set = set(protein_ids)
        with open(files[0], 'r') as f:
            for line in f:
                protein_id = line.split('\t', 1)[0]
                if protein_id in protein_id_set:
                    _, *preds = line.strip().split('\t')
                    predictions[protein_id] = {
                        go_id: float(score_str)
                        for pred in preds
                        for go_id, score_str in [pred.split('|', 1)]
                    }
    return predictions


def load_predictions_from_file(
    predictions_file: str, test_proteins: Dict[str, List[str]]
) -> Dict[str, Dict[str, float]]:
    """Load predictions from single TSV file."""
    protein_id_set = set()
    for ids in test_proteins.values():
        protein_id_set.update(ids)
    predictions = {}
    if not os.path.exists(predictions_file):
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


def load_optimized_from_dir(
    optimized_dir: str, test_proteins: Dict[str, List[str]]
) -> Dict[str, Dict[str, float]]:
    """Load optimized annotations from per-taxon TSV files (optimized_*_taxon_*.tsv).

    Original scores are preserved so that the optimized set can be evaluated at the
    same Fmax threshold as the raw predictions, making the before/after comparison
    apples-to-apples.
    """
    optimized = defaultdict(dict)
    for taxon_id, protein_ids in test_proteins.items():
        files = glob.glob(os.path.join(optimized_dir, f'optimized_*_taxon_{taxon_id}.tsv'))
        if not files:
            continue
        protein_id_set = set(protein_ids)
        with open(files[0], 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if not parts:
                    continue
                protein_id = parts[0]
                if protein_id in protein_id_set:
                    for pred in parts[1:]:
                        if '|' in pred:
                            go_id, score_str = pred.split('|', 1)
                            score = float(score_str)
                            if score > 0:
                                optimized[protein_id][go_id] = score
    return dict(optimized)


def load_optimized_from_file(
    optimized_file: str, test_proteins: Dict[str, List[str]]
) -> Dict[str, Dict[str, float]]:
    """Load optimized annotations from single TSV file.

    Original scores are preserved so that the optimized set can be evaluated at the
    same Fmax threshold as the raw predictions.
    """
    protein_id_set = set()
    for ids in test_proteins.values():
        protein_id_set.update(ids)
    optimized = {}
    if not os.path.exists(optimized_file):
        return optimized
    with open(optimized_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if not parts:
                continue
            protein_id = parts[0]
            if protein_id in protein_id_set:
                optimized[protein_id] = {}
                for pred in parts[1:]:
                    if '|' in pred:
                        go_id, score_str = pred.split('|', 1)
                        score = float(score_str)
                        if score > 0:
                            optimized[protein_id][go_id] = score
    return optimized


def parse_bin_spec(bin_spec: str) -> List[tuple]:
    """Parse bin specification string into list of (low, high, label) tuples.

    Examples: "1" -> (1,1), "2-5" -> (2,5), "101+" -> (101, inf)
    """
    bins = []
    for token in bin_spec.split(','):
        token = token.strip()
        if token.endswith('+'):
            low = int(token[:-1])
            bins.append((low, float('inf'), token))
        elif '-' in token:
            low, high = token.split('-', 1)
            bins.append((int(low), int(high), token))
        else:
            val = int(token)
            bins.append((val, val, token))
    return bins


def compute_binned_summary(rows: List[dict], bin_spec: str = "1,2-5,6-10,11-25,26-50,51-100,101+") -> List[dict]:
    """Group terms by annotation count bins and compute summary stats per bin."""
    bins = parse_bin_spec(bin_spec)
    bin_rows = {label: [] for _, _, label in bins}

    for r in rows:
        n = r['n_annotations']
        for low, high, label in bins:
            if low <= n <= high:
                bin_rows[label].append(r)
                break

    summary = []
    for low, high, label in bins:
        group = bin_rows[label]
        n_terms = len(group)
        if n_terms == 0:
            summary.append({'bin': label, 'n_terms': 0})
            continue

        auc_vals = [r['auc_pred'] for r in group if not math.isnan(r['auc_pred'])]
        ba_vals = [r['balanced_acc_opt'] for r in group]
        diff_sens = [r['diff_sensitivity'] for r in group]
        diff_spec = [r['diff_specificity'] for r in group]

        summary.append({
            'bin': label,
            'n_terms': n_terms,
            'mean_auc_pred': statistics.mean(auc_vals) if auc_vals else float('nan'),
            'median_auc_pred': statistics.median(auc_vals) if auc_vals else float('nan'),
            'std_auc_pred': statistics.stdev(auc_vals) if len(auc_vals) > 1 else 0.0,
            'mean_balanced_acc_opt': statistics.mean(ba_vals),
            'mean_diff_sensitivity': statistics.mean(diff_sens),
            'mean_diff_specificity': statistics.mean(diff_spec),
            'pct_sens_improved': sum(1 for d in diff_sens if d > 0) / n_terms * 100,
            'pct_spec_improved': sum(1 for d in diff_spec if d > 0) / n_terms * 100,
        })

    return summary


def plot_auc_by_annotation_count(binned_summary: List[dict], output_path: str, title: str = ""):
    """Generate a grouped bar chart of mean AUC (pred) vs balanced accuracy (opt) per bin."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    labels = [b['bin'] for b in binned_summary if b['n_terms'] > 0]
    pred_vals = [b['mean_auc_pred'] for b in binned_summary if b['n_terms'] > 0]
    opt_vals = [b['mean_balanced_acc_opt'] for b in binned_summary if b['n_terms'] > 0]

    x = range(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar([i - width / 2 for i in x], pred_vals, width, label='Predictions (AUC)')
    bars2 = ax.bar([i + width / 2 for i in x], opt_vals, width, label='Optimized (Balanced Acc.)')

    ax.set_xlabel('Number of Annotations')
    ax.set_ylabel('Performance')
    ax.set_title(title or 'Term-Centric Performance by Annotation Count')
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.legend()
    ax.set_ylim(0, 1.05)

    # Add term counts above bars
    for i, b in enumerate([b for b in binned_summary if b['n_terms'] > 0]):
        ax.text(i, max(pred_vals[i], opt_vals[i]) + 0.02, f"n={b['n_terms']}",
                ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Saved plot to {output_path}")


def run_report(
    predictions: Dict[str, Dict[str, float]],
    optimized: Dict[str, Dict[str, float]],
    test_proteins: Dict[str, List[str]],
    annotations_dir: str,
    go_file: str,
    subontology: str,
    output_file: str,
    top_k: int,
    output_dir: Optional[str] = None,
    min_annotations: int = 0,
    bin_spec: str = "1,2-5,6-10,11-25,26-50,51-100,101+",
    plot: bool = False,
) -> None:
    """Compute per-term metrics, rank by difference, write report."""
    if PIPELINE_DIR not in sys.path:
        sys.path.append(PIPELINE_DIR)
    if BASE_DIR not in sys.path:
        sys.path.append(BASE_DIR)
    if GAEF_DIR not in sys.path:
        sys.path.append(GAEF_DIR)

    import importlib.util
    genome_scores_path = os.path.join(PIPELINE_DIR, 'genome_scores_evaluator.py')
    spec = importlib.util.spec_from_file_location("genome_scores_evaluator", genome_scores_path)
    genome_scores = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(genome_scores)

    load_annotations = genome_scores.load_annotations
    Ontology = genome_scores.Ontology
    compute_per_term_metrics = genome_scores.compute_per_term_metrics
    compute_per_term_auc = genome_scores.compute_per_term_auc
    propagate_predictions_dict = genome_scores.propagate_predictions_dict

    if not os.path.exists(annotations_dir):
        print(f"[Error] Annotations directory not found: {annotations_dir}")
        return
    if not os.path.exists(go_file):
        print(f"[Error] GO file not found: {go_file}")
        return

    protein_to_taxon = {}
    for taxon_id, protein_ids in test_proteins.items():
        for protein_id in protein_ids:
            protein_to_taxon[protein_id] = taxon_id

    pooled_predictions = {}
    for protein_id, scores in predictions.items():
        taxon_id = protein_to_taxon.get(protein_id)
        if taxon_id is not None:
            pooled_predictions[f"{taxon_id}:{protein_id}"] = scores

    pooled_optimized = {}
    for protein_id, terms_dict in optimized.items():
        taxon_id = protein_to_taxon.get(protein_id)
        if taxon_id is not None:
            pooled_optimized[f"{taxon_id}:{protein_id}"] = terms_dict

    pooled_annotations = {}
    for taxon_id, protein_ids in test_proteins.items():
        annot_file = os.path.join(annotations_dir, f'annots_taxon_{taxon_id}.tsv')
        if not os.path.exists(annot_file):
            continue
        taxon_annotations = load_annotations(annot_file)
        protein_set = set(protein_ids)
        for prot_id, go_terms in taxon_annotations.items():
            if prot_id in protein_set:
                pooled_annotations[f"{taxon_id}:{prot_id}"] = go_terms

    if not pooled_annotations:
        print("[Error] No annotations found for provided test proteins")
        return

    go = Ontology(go_file)

    pooled_predictions_prop = propagate_predictions_dict(pooled_predictions, go, subontology)
    pooled_optimized_prop = propagate_predictions_dict(pooled_optimized, go, subontology)

    per_term_pred, tmax = compute_per_term_metrics(
        go, subontology, pooled_annotations, pooled_predictions_prop, threshold=None
    )
    if per_term_pred is None:
        print("[Error] Could not compute per-term metrics for predictions")
        return

    # Evaluate optimized at the same Fmax threshold so the comparison is apples-to-apples.
    # The adjustment is removal-only (original scores preserved), so any term present in
    # optimized but absent in predictions was below tmax and stays below tmax here.
    per_term_opt, _ = compute_per_term_metrics(
        go, subontology, pooled_annotations, pooled_optimized_prop, threshold=tmax
    )
    if per_term_opt is None:
        print("[Error] Could not compute per-term metrics for optimized")
        return

    # Per-term AUC for predictions (continuous scores)
    per_term_auc_pred = compute_per_term_auc(
        go, subontology, pooled_annotations, pooled_predictions_prop
    )

    total_proteins = len(pooled_annotations)

    rows = []
    for term_id in per_term_pred:
        p = per_term_pred[term_id]
        o = per_term_opt.get(term_id, {'sensitivity': 0.0, 'specificity': 0.0})
        # Signed differences: positive means solver improved
        diff_sens = o['sensitivity'] - p['sensitivity']
        diff_spec = o['specificity'] - p['specificity']
        rank_score = max(abs(diff_sens), abs(diff_spec))

        n_ann = p['n_annotations']
        n_neg = total_proteins - n_ann

        tp_pred = round(p['sensitivity'] * n_ann)
        tp_opt = round(o['sensitivity'] * n_ann)
        fp_pred = round((1 - p['specificity']) * n_neg) if n_neg > 0 else 0
        fp_opt = round((1 - o['specificity']) * n_neg) if n_neg > 0 else 0

        delta_tp = tp_opt - tp_pred
        delta_fp = fp_opt - fp_pred

        # Per-term AUC (predictions) and balanced accuracy (optimized)
        auc_pred = float('nan')
        if per_term_auc_pred and term_id in per_term_auc_pred:
            auc_pred = per_term_auc_pred[term_id]['auc']
        balanced_acc_opt = (o['sensitivity'] + o['specificity']) / 2.0

        rows.append({
            'term_id': term_id,
            'n_annotations': n_ann,
            'auc_pred': auc_pred,
            'balanced_acc_opt': balanced_acc_opt,
            'sensitivity_pred': p['sensitivity'],
            'sensitivity_opt': o['sensitivity'],
            'diff_sensitivity': diff_sens,
            'specificity_pred': p['specificity'],
            'specificity_opt': o['specificity'],
            'diff_specificity': diff_spec,
            'rank_score': rank_score,
            'tp_pred': tp_pred,
            'tp_opt': tp_opt,
            'delta_tp': delta_tp,
            'fp_pred': fp_pred,
            'fp_opt': fp_opt,
            'delta_fp': delta_fp,
        })

    # Flip-precision summary (printed before any filtering).
    # Since optimized is evaluated at the same tmax as predictions, the adjustment is
    # removal-only: sensitivity can only decrease (TP lost) and specificity can only
    # increase (FP removed). Flip precision = FPs removed / total terms changed.
    changed = [r for r in rows if r['rank_score'] > 0]
    tps_lost  = sum(1 for r in changed if r['delta_tp'] < 0)
    fps_removed = sum(1 for r in changed if r['delta_fp'] < 0)
    total_changed = tps_lost + fps_removed
    flip_precision = fps_removed / total_changed if total_changed > 0 else float('nan')
    unchanged = sum(1 for r in rows if r['rank_score'] == 0)
    print(f"Terms changed by solver: {total_changed}  "
          f"(FPs removed: {fps_removed}, TPs lost: {tps_lost})  "
          f"Flip precision: {flip_precision:.3f}  "
          f"Terms unchanged: {unchanged}")
    print(f"F-max threshold (tmax): {tmax}")

    # Apply min_annotations filter
    if min_annotations > 0:
        rows = [r for r in rows if r['n_annotations'] >= min_annotations]

    rows.sort(key=lambda r: r['rank_score'], reverse=True)
    rows_with_rank_score_greater_than_0 = [r for r in rows if r['rank_score'] > 0]
    top_rows = rows_with_rank_score_greater_than_0[:top_k]

    # Write top-K report (backward-compatible output)
    os.makedirs(os.path.dirname(os.path.abspath(output_file)) or '.', exist_ok=True)
    with open(output_file, 'w') as f:
        cols = ['term_id', 'n_annotations', 'auc_pred', 'balanced_acc_opt',
                'sensitivity_pred', 'sensitivity_opt', 'diff_sensitivity',
                'specificity_pred', 'specificity_opt', 'diff_specificity', 'rank_score',
                'tp_pred', 'tp_opt', 'delta_tp', 'fp_pred', 'fp_opt', 'delta_fp']
        f.write('\t'.join(cols) + '\n')
        for r in top_rows:
            f.write('\t'.join(str(r[c]) for c in cols) + '\n')

    print(f"Total test proteins: {total_proteins}")
    print(f"F-max threshold for predictions: {tmax}")
    print(f"Wrote top {len(top_rows)} terms to {output_file}")

    # Extended outputs when output_dir is specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        # Full per-term detail TSV
        detail_file = os.path.join(output_dir, f'per_term_detail_{subontology}.tsv')
        detail_cols = ['term_id', 'n_annotations', 'auc_pred', 'balanced_acc_opt',
                       'sensitivity_pred', 'sensitivity_opt', 'diff_sensitivity',
                       'specificity_pred', 'specificity_opt', 'diff_specificity', 'rank_score',
                       'tp_pred', 'tp_opt', 'delta_tp', 'fp_pred', 'fp_opt', 'delta_fp']
        with open(detail_file, 'w') as f:
            f.write('\t'.join(detail_cols) + '\n')
            for r in rows:
                f.write('\t'.join(str(r[c]) for c in detail_cols) + '\n')
        print(f"Wrote {len(rows)} terms (full detail) to {detail_file}")

        # Binned summary TSV
        binned = compute_binned_summary(rows, bin_spec)
        summary_file = os.path.join(output_dir, f'binned_summary_{subontology}.tsv')
        summary_cols = ['bin', 'n_terms', 'mean_auc_pred', 'median_auc_pred', 'std_auc_pred',
                        'mean_balanced_acc_opt', 'mean_diff_sensitivity', 'mean_diff_specificity',
                        'pct_sens_improved', 'pct_spec_improved']
        with open(summary_file, 'w') as f:
            f.write('\t'.join(summary_cols) + '\n')
            for b in binned:
                vals = [str(b.get(c, '')) for c in summary_cols]
                f.write('\t'.join(vals) + '\n')
        print(f"Wrote binned summary to {summary_file}")

        # Optional plot
        if plot:
            plot_file = os.path.join(output_dir, f'auc_by_annotation_count_{subontology}.png')
            plot_auc_by_annotation_count(
                binned, plot_file,
                title=f'Term-Centric Performance by Annotation Count ({subontology.upper()})'
            )


def main(
    predictions_dir: Optional[str],
    predictions_file: Optional[str],
    optimized_dir: Optional[str],
    optimized_file: Optional[str],
    test_proteins_file: str,
    output_file: str,
    annotations_dir: str,
    go_file: str,
    subontology: str,
    top_k: int,
    output_dir: Optional[str] = None,
    min_annotations: int = 0,
    bin_spec: str = "1,2-5,6-10,11-25,26-50,51-100,101+",
    plot: bool = False,
) -> None:
    if predictions_file and predictions_dir:
        raise ValueError("Provide either --predictions_file or --predictions_dir, not both")
    if not predictions_file and not predictions_dir:
        raise ValueError("Provide either --predictions_file or --predictions_dir")
    if optimized_file and optimized_dir:
        raise ValueError("Provide either --optimized_file or --optimized_dir, not both")
    if not optimized_file and not optimized_dir:
        raise ValueError("Provide either --optimized_file or --optimized_dir")

    test_proteins = load_test_proteins(test_proteins_file)

    if predictions_file:
        predictions = load_predictions_from_file(predictions_file, test_proteins)
        print(f"Loaded predictions from {predictions_file}")
    else:
        predictions = load_predictions_from_dir(predictions_dir, test_proteins)
        print(f"Loaded predictions from {predictions_dir}")

    if optimized_file:
        optimized = load_optimized_from_file(optimized_file, test_proteins)
        print(f"Loaded optimized from {optimized_file}")
    else:
        optimized = load_optimized_from_dir(optimized_dir, test_proteins)
        print(f"Loaded optimized from {optimized_dir}")

    run_report(
        predictions=predictions,
        optimized=optimized,
        test_proteins=test_proteins,
        annotations_dir=annotations_dir,
        go_file=go_file,
        subontology=subontology,
        output_file=output_file,
        top_k=top_k,
        output_dir=output_dir,
        min_annotations=min_annotations,
        bin_spec=bin_spec,
        plot=plot,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Report top terms by metric difference (predictions vs optimized annotations)'
    )
    parser.add_argument('--predictions_dir', type=str, default=None)
    parser.add_argument('--predictions_file', type=str, default=None)
    parser.add_argument('--optimized_dir', type=str, default=None)
    parser.add_argument('--optimized_file', type=str, default=None)
    parser.add_argument('--test_proteins_file', type=str, default='proteins_by_date_23-MAY-2024.tsv')
    parser.add_argument('--output_file', type=str, default='term_metric_diff_report.tsv')
    parser.add_argument('--annotations_dir', type=str,
                        default="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic")
    parser.add_argument('--go_file', type=str, default="data/go-basic.obo")
    parser.add_argument('--subontology', type=str, default='cc', choices=['cc', 'bp', 'mf'])
    parser.add_argument('--top_k', type=int, default=50)
    parser.add_argument('--output_dir', type=str, default=None,
                        help='Directory for extended outputs (full detail TSV, binned summary, plot)')
    parser.add_argument('--min_annotations', type=int, default=0,
                        help='Minimum annotation count to include a term (default: 0, no filter)')
    parser.add_argument('--bins', type=str, default="1,2-5,6-10,11-25,26-50,51-100,101+",
                        help='Annotation count bin specification (default: 1,2-5,6-10,11-25,26-50,51-100,101+)')
    parser.add_argument('--plot', action='store_true',
                        help='Generate matplotlib plot (requires --output_dir)')
    args = parser.parse_args()

    main(
        args.predictions_dir,
        args.predictions_file,
        args.optimized_dir,
        args.optimized_file,
        args.test_proteins_file,
        args.output_file,
        args.annotations_dir,
        args.go_file,
        args.subontology,
        args.top_k,
        output_dir=args.output_dir,
        min_annotations=args.min_annotations,
        bin_spec=args.bins,
        plot=args.plot,
    )
