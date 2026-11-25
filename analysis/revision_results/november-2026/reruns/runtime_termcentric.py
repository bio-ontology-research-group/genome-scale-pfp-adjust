"""Runtimes and term-centric flip ranking from the accepted corrected task outputs.

Reads out/rev5-corrected-baselines/<predictor>/<aspect>/<taxon>.json.gz written by
analysis/revision_audit_run.py. Stage 1 seconds (model build, solve and export)
are recorded per run as seconds_stage1; Stage 2 solve seconds are in
runs[*].solver.seconds. Term-centric flips compare 'before' and 'stage1'
test_scores at the native threshold on the eligible test proteins, restricted
to pairs whose original score exceeds 0.1 (as in the original analysis).
"""
import argparse
import gzip
import json
import statistics as st
import sys
from collections import defaultdict, Counter
from pathlib import Path

PRED = ['seq_sim_timeset', 'deepgo-se', 'sprof-go']
NATIVE = {'deepgo-se': .3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--code', required=True)
    ap.add_argument('--go', required=True)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology
    go = Ontology(a.go)
    def name(g):
        try:
            return go.get_term(g)['name']
        except Exception:
            return ''
    res = {'runtime': {}, 'term_centric': {}}
    flipped_by = defaultdict(set)
    for pred in PRED:
        tau = NATIVE.get(pred, .001)
        for asp in ['cc', 'mf', 'bp']:
            s1, s1u, s2 = [], [], defaultdict(list)
            cost, count = Counter(), Counter()
            n_flips = 0
            for f in sorted(Path(a.out_dir, pred, asp).glob('*.json.gz')):
                d = json.load(gzip.open(f, 'rt'))
                runs = [r for r in d['runs'] if abs(r['tau'] - tau) < 1e-12]
                by = {r['mode']: r for r in runs}
                s1.append(by['stage1']['seconds_stage1'])
                for m in ['stage2_standard', 'stage2_demotion']:
                    if m in by and by[m].get('solver'):
                        s2[m].append((by[m]['solver']['seconds'], by[m]['solver']['status']))
                b, s = by['before']['test_scores'], by['stage1']['test_scores']
                for p, scores in b.items():
                    after = s.get(p, {})
                    for g, v in scores.items():
                        if v > .1 and after.get(g, 0.0) != v:
                            cost[g] += abs(v - after.get(g, 0.0))
                            count[g] += 1
                            n_flips += 1
            key = f'{pred}_{asp}'
            res['runtime'][key] = {
                'n': len(s1), 'stage1_mean': st.mean(s1), 'stage1_sd': st.pstdev(s1), 'stage1_max': max(s1),
                **{f'{m}_mean': st.mean(x for x, _ in v) for m, v in s2.items()},
                **{f'{m}_sd': st.pstdev([x for x, _ in v]) for m, v in s2.items()},
                **{f'{m}_max': max(x for x, _ in v) for m, v in s2.items()},
                **{f'{m}_status': dict(Counter(y for _, y in v)) for m, v in s2.items()}}
            total = sum(cost.values())
            top = cost.most_common(10)
            for g in cost:
                flipped_by[(asp, g)].add(pred)
            res['term_centric'][key] = {
                'n_flips': n_flips, 'total_cost': total,
                'top10_share': sum(c for _, c in top) / total if total else None,
                'top10': [{'term': g, 'name': name(g), 'cost': c, 'flips': count[g]} for g, c in top]}
    res['flipped_by_all_three'] = {asp: sorted(f'{g} {name(g)}' for (x, g), ps in flipped_by.items()
                                              if x == asp and len(ps) == 3) for asp in ['cc', 'mf', 'bp']}
    Path(a.output).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1)[:6000])


if __name__ == '__main__':
    main()
