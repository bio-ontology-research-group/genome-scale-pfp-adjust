"""Per-term score summaries of a proteome's canonical, propagated predictions.

For each GO term with any protein score > tau, record the number of proteins and
the score sum in fixed score bins, plus the maximum score. Any thresholded
summed-margin cost sum_i (s_ij - t)_+ for t in BIN_EDGES, the max score and
protein counts can be reconstructed exactly from this summary.
"""
import argparse
import gzip
import json
import sys
from pathlib import Path

BIN_EDGES = [0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--resources', required=True)
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--aspect', required=True)
    ap.add_argument('--tau', type=float, required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology, NAMESPACES
    from taxon_consistency import adjust_ortools as tc
    go = Ontology(Path(a.resources) / 'go-basic.obo')
    raw = go.propagate_scores(go.normalize_scores(tc.load_predictions(a.predictions), NAMESPACES[a.aspect], {}), {})
    nb = len(BIN_EDGES)
    stats = {}
    for p, sc in raw.items():
        for g, s in sc.items():
            if s <= a.tau:
                continue
            st = stats.get(g)
            if st is None:
                st = stats[g] = {'n': [0] * nb, 's': [0.0] * nb, 'max': 0.0}
            b = 0
            while b + 1 < nb and s >= BIN_EDGES[b + 1]:
                b += 1
            st['n'][b] += 1
            st['s'][b] += s
            if s > st['max']:
                st['max'] = s
    margin = {g: sum(st['s']) - a.tau * sum(st['n']) for g, st in stats.items()}
    with gzip.open(a.output, 'wt') as fh:
        json.dump({'tau': a.tau, 'aspect': a.aspect, 'predictions': a.predictions, 'n_proteins': len(raw),
                   'bin_edges': BIN_EDGES, 'terms': stats, 'margin': margin}, fh)
    print('terms', len(stats), 'proteins', len(raw))


if __name__ == '__main__':
    main()
