"""Summarize lineage_variants.py outputs (and the accepted out-r1 robust records)."""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

RANKS = ['none', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def fmt(c):
    return ' '.join(f'{r[:3]}{c.get(r, 0)}' for r in RANKS if c.get(r, 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--baseline', help='extensions out-r1/lineage dir (robust_evaluation)')
    ap.add_argument('--ceiling')
    ap.add_argument('--output')
    a = ap.parse_args()
    groups = defaultdict(list)
    for f in sorted(Path(a.dir).glob('*.json')):
        pred = f.stem.rsplit('_', 1)[0]
        for r in json.loads(f.read_text())['records']:
            groups[(r['variant'], pred, r['aspect'])].append(r)
    if a.baseline:
        for pdir in Path(a.baseline).iterdir():
            for f in pdir.glob('*.json'):
                d = json.loads(f.read_text())
                for rec in d['records']:
                    ev = rec['robust_evaluation']
                    cat = ('unresolved' if not ev['inferred_ranks'] else
                           'correct' if ev['is_exact_prefix'] else 'wrong')
                    groups[('outr1_robust', pdir.name, rec['aspect'])].append(
                        dict(taxon=d['taxon'], category=cat, deepest=ev['deepest_correct_rank'],
                             divergence=ev['first_divergence_rank']))
    ceil = {}
    if a.ceiling:
        ceil = {r['taxon']: r for r in json.loads(Path(a.ceiling).read_text())['rows']}
    table = []
    for key in sorted(groups):
        rs = groups[key]
        cats = Counter(r['category'] for r in rs)
        deep = Counter(r['deepest'] for r in rs if r['category'] == 'correct')
        div = Counter(r['divergence'] for r in rs if r['category'] == 'wrong')
        row = dict(variant=key[0], predictor=key[1], aspect=key[2], n=len(rs),
                   correct=cats['correct'], wrong=cats['wrong'], unresolved=cats['unresolved'],
                   other=len(rs) - cats['correct'] - cats['wrong'] - cats['unresolved'],
                   deepest_correct=dict(deep), divergence_wrong=dict(div))
        if ceil:
            wrong_dom = Counter(ceil[r['taxon']]['domain'] for r in rs if r['category'] == 'wrong'
                                and r['taxon'] in ceil)
            row['wrong_by_true_domain'] = dict(wrong_dom)
            # correct at the ceiling (deepest achievable rank)
            row['at_ceiling'] = sum(1 for r in rs if r['category'] == 'correct' and r['taxon'] in ceil
                                    and r['deepest'] == ceil[r['taxon']]['ceiling'])
        table.append(row)
        print(f"{key[0]:14s} {key[1]:16s} {key[2]:5s} n={len(rs):3d} corr={cats['correct']:3d} "
              f"wrong={cats['wrong']:3d} unres={cats['unresolved']:3d} oth={row['other']:2d} "
              f"| deep: {fmt(deep):40s} | div: {fmt(div)} | {row.get('wrong_by_true_domain', '')} "
              f"atceil={row.get('at_ceiling', '')}")
    if a.output:
        Path(a.output).write_text(json.dumps(table, indent=1))


if __name__ == '__main__':
    main()
