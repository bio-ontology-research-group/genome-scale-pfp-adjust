"""Characterise how the ORIGINAL pipeline's optimized/ files differ from predictions/.

For each organism: number of raw (protein, term) entries, entries removed, removed entries
with score > tau and <= tau, number of distinct terms removed, and whether removed terms were
removed at all scores (term-level) or only for above-tau proteins. Also counts active pairs
(> tau) as the original flip table would.
"""
import argparse
import glob
import json
import os
import re
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument('--root', required=True)
ap.add_argument('--pred', default='deepgo-se')
ap.add_argument('--aspect', default='cc')
ap.add_argument('--tau', type=float, default=0.3)
ap.add_argument('--opt', default='optimized')
ap.add_argument('--output', required=True)
a = ap.parse_args()


def load(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0]:
                continue
            d = {}
            for tok in parts[1:]:
                if '|' in tok:
                    g, s = tok.split('|', 1)
                    d[g] = float(s)
            out[parts[0]] = d
    return out

rows = []
tot = Counter()
base = f'{a.root}/{a.pred}_results/{a.aspect}'
for pf in sorted(glob.glob(f'{base}/predictions/predictions_fold_*_taxon_*.tsv')):
    taxon = re.search(r'taxon_(\d+)', pf).group(1)
    of = pf.replace('/predictions/', f'/{a.opt}/').replace('predictions_fold', 'optimized_fold')
    if not os.path.exists(of):
        continue
    P, O = load(pf), load(of)
    c = Counter()
    removed_terms = Counter(); kept_terms_with_removed = set()
    changed_score = 0
    for p, sc in P.items():
        o = O.get(p, {})
        for g, s in sc.items():
            c['raw'] += 1
            if s > a.tau:
                c['active'] += 1
            if g not in o:
                c['removed'] += 1
                c['removed_above' if s > a.tau else 'removed_below'] += 1
                removed_terms[g] += 1
            elif o[g] != s and abs(o[g] - s) > 1e-6:
                changed_score += 1
        for g in o:
            if g not in sc:
                c['added'] += 1
    # term-level: removed terms still present anywhere in optimized?
    present_after = set()
    for p, o in O.items():
        present_after.update(o)
    partial = sorted(g for g in removed_terms if g in present_after)
    c['changed_score'] = changed_score
    c['n_removed_terms'] = len(removed_terms)
    c['n_removed_terms_still_present'] = len(partial)
    c['proteins_touched'] = sum(1 for p, sc in P.items() if any(g not in O.get(p, {}) for g in sc))
    rows.append({'taxon': taxon, **c, 'top_removed': removed_terms.most_common(8)})
    tot.update({k: v for k, v in c.items()})
    print(taxon, dict(c), flush=True)
json.dump({'total': tot, 'rows': rows}, open(a.output, 'w'), indent=1)
print('TOTAL', dict(tot))
