"""Diagnose wrong lineages: optimum vs true-lineage repair under a variant."""
import argparse
import json
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

HIST = ['<.01', '.01-.1', '.1-.3', '.3-.5', '.5-.7', '.7-.9', '>=.9']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', nargs='+', required=True)
    ap.add_argument('--ceiling', required=True)
    ap.add_argument('--variant', default='margin')
    ap.add_argument('--aspect', default='joint')
    a = ap.parse_args()
    ceil = {r['taxon']: r for r in json.loads(Path(a.ceiling).read_text())['rows']}
    by = defaultdict(list)
    for d in a.dir:
        for f in Path(d).glob('*.json'):
            pred = f.stem.rsplit('_', 1)[0]
            for r in json.loads(f.read_text())['records']:
                if r['variant'] == a.variant and r['aspect'] == a.aspect and 'truth_cost' in r:
                    by[pred].append(r)
    for pred, rs in sorted(by.items()):
        wrong = [r for r in rs if r['category'] == 'wrong']
        print(f'== {pred} {a.variant} {a.aspect}: n={len(rs)} wrong={len(wrong)}')
        for dom in ['Bacteria', 'Archaea', 'Eukaryota']:
            w = [r for r in wrong if ceil[r['taxon']]['domain'] == dom]
            if not w:
                continue
            ratio = [r['truth_cost'] / max(r['demoted_cost'], 1e-9) for r in w]
            tn = [r['truth_n_demoted'] for r in w]
            on = [r['n_demoted'] for r in w]
            th = [sum(x) for x in zip(*[r['truth_demoted_max_hist'] for r in w])]
            oh = [sum(x) for x in zip(*[r['demoted_max_hist'] for r in w])]
            top = Counter(g for r in w for g, *_ in r['truth_top_demoted'][:5])
            otop = Counter(g for r in w for g, *_ in r['top_demoted'][:5])
            print(f'  {dom}: wrong={len(w)} truth/opt cost median={st.median(ratio):.2f} '
                  f'(range {min(ratio):.2f}-{max(ratio):.2f}); demoted terms truth median={st.median(tn)} '
                  f'opt median={st.median(on)}')
            print('     truth-demoted max-score hist', dict(zip(HIST, th)))
            print('     opt-demoted   max-score hist', dict(zip(HIST, oh)))
            print('     top truth-demoted', top.most_common(8))
            print('     top opt-demoted  ', otop.most_common(8))


if __name__ == '__main__':
    main()
