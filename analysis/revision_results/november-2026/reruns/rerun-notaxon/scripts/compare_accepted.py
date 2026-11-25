"""Compare rerun outputs to accepted rev5 outputs per organism (before/stage1/stage1_unknown at native tau)."""
import gzip, json, sys
from pathlib import Path
new_root, acc_root = Path(sys.argv[1]), Path(sys.argv[2])
tasks = [l.split() for l in Path(sys.argv[3]).read_text().splitlines() if l.strip()]
keys = ['n_proteins', 'active_predictions', 'n_active_terms', 'demotions', 'promotions', 'known_taxon_consistent',
        'gaef_taxon_consistent', 'violating_terms', 'process_coherence', 'pathway_coherence', 'essential_completeness', 'complex_coherence']
bad = 0
for pred, asp, t in tasks:
    a = acc_root / pred / asp / f'{t}.json.gz'; n = new_root / pred / asp / f'{t}.json.gz'
    if not (a.exists() and n.exists()):
        print('MISSING', pred, asp, t, a.exists(), n.exists()); bad += 1; continue
    A = json.load(gzip.open(a, 'rt')); N = json.load(gzip.open(n, 'rt'))
    tau = A['native_tau']
    msgs = []
    for mode in ['before', 'stage1', 'stage1_unknown']:
        ra = [r for r in A['runs'] if r['mode'] == mode and r['tau'] == tau]
        rn = [r for r in N['runs'] if r['mode'] == mode and r['tau'] == tau]
        ra, rn = ra[0], rn[0]
        for k in keys:
            if ra.get(k) != rn.get(k):
                msgs.append(f'{mode}.{k}: {str(ra.get(k))[:60]} vs {str(rn.get(k))[:60]}')
        if mode == 'stage1_unknown':
            if sorted(ra['solver_taxons']) != sorted(rn['solver_taxons']):
                msgs.append('solver_taxons differ')
            if ra['retained_active_terms'] != rn['retained_active_terms']:
                msgs.append('retained terms differ')
        common = set(ra['test_scores']) & set(rn['test_scores'])
        diff = sum(ra['test_scores'][p] != rn['test_scores'][p] for p in common)
        if diff:
            msgs.append(f'{mode}: test_scores differ on {diff}/{len(common)}')
        if mode != 'before':
            ma = set(ra.get('modified_test_proteins', [])) & common; mn = set(rn.get('modified_test_proteins', [])) & common
            if ma != mn:
                msgs.append(f'{mode}: modified proteins {len(ma)} vs {len(mn)}')
    print('OK  ' if not msgs else 'DIFF', pred, asp, t, '; '.join(msgs))
    bad += bool(msgs)
print('n_bad', bad, 'of', len(tasks))
