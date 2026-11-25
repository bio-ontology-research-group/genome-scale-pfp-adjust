"""Reproduce the ORIGINAL-submission micro-averaged CAFA metrics with the original
evaluator code ($BASE/code == commit 3ffb56a pipeline/) on the original inputs
(/ibex/scratch/projects/c2014/asaad/swissprot_proteomes_folds/<pred>_results/<asp>/{predictions,optimized}),
and capture the PR points that the original AUPR (np.trapz over argsort(recall)) integrated.

For each condition, writes the 101 (threshold, precision, recall) points and AUPR computed
(i) as the original code did, (ii) with a stable sort, (iii) with the corrected rule
(empty-prediction endpoint (0,1), max precision per repeated recall).
"""
import argparse
import json
import sys
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--oldcode', required=True)
ap.add_argument('--gaef', required=True)
ap.add_argument('--go', required=True)
ap.add_argument('--root', required=True)
ap.add_argument('--pred', default='deepgo-se')
ap.add_argument('--aspect', default='cc')
ap.add_argument('--modes', default='predictions,optimized')
ap.add_argument('--taxa', default=None, help='optional file listing taxa to restrict pooling')
ap.add_argument('--output', required=True)
a = ap.parse_args()
sys.path[:0] = [a.oldcode, a.oldcode + '/pipeline', a.gaef]
import evaluate_directory as ed  # old code
import genome_scores_evaluator as gse

captured = []
_orig_trapz = np.trapz


def spy(y, x=None, *args, **kw):
    captured.append((np.array(y, dtype=float).copy(), np.array(x, dtype=float).copy()))
    return _orig_trapz(y, x, *args, **kw)


np.trapz = spy
gse.np.trapz = spy

# Capture unsorted (threshold-ordered) precision/recall too by wrapping the
# per-threshold evaluator.
order = []
_orig_eval = gse.evaluate_annotations_optimized


def eval_spy(go, labels, preds, ic_cache=None):
    res = _orig_eval(go, labels, preds, ic_cache)
    order.append((res[1], res[2], sum(1 for p in preds if p)))
    return res


gse.evaluate_annotations_optimized = eval_spy

if a.taxa:
    keep = set(open(a.taxa).read().split())
    _orig_find = ed.find_prediction_files

    def find(d, pattern='*.tsv'):
        out = []
        for f in _orig_find(d, pattern):
            ident = ed.extract_identifier_from_filename(f.split('/')[-1])
            t = ident[1] if isinstance(ident, tuple) else ident
            if t in keep:
                out.append(f)
        return out
    ed.find_prediction_files = find


def aupr_new(prec, rec):
    curve = {0.0: 1.0}
    for r, p in zip(rec, prec):
        curve[r] = max(curve.get(r, 0.0), p)
    pts = sorted(curve.items())
    return sum((r1 - r0) * (p0 + p1) / 2 for (r0, p0), (r1, p1) in zip(pts, pts[1:]))


result = {}
for mode in a.modes.split(','):
    captured.clear(); order.clear()
    d = f'{a.root}/{a.pred}_results/{a.aspect}/{mode}'
    m = ed.compute_micro_averaged_metrics(d, f'{a.root}/annotations-go-basic', a.go, a.aspect)
    pts = order[-101:]
    prec = [p for p, r, n in pts]
    rec = [r for p, r, n in pts]
    nprot = [n for p, r, n in pts]
    P, R = np.array(prec), np.array(rec)
    idx = np.argsort(R, kind='stable')
    result[mode] = {
        'metrics': m,
        'thresholds': [t / 100 for t in range(101)],
        'precision': prec, 'recall': rec, 'n_proteins_with_pred': nprot,
        'aupr_original_code': float(_orig_trapz(captured[-1][0], captured[-1][1])) if captured else None,
        'aupr_stable_sort': float(_orig_trapz(P[idx], R[idx])),
        'aupr_threshold_order': float(-_orig_trapz(P, R)),
        'aupr_corrected_rule': aupr_new(prec, rec),
        'aupr_corrected_rule_excl_t0': aupr_new(prec[1:], rec[1:]),
    }
    print(mode, {k: v for k, v in result[mode].items() if k.startswith('aupr')}, flush=True)
json.dump(result, open(a.output, 'w'), indent=1, default=float)
