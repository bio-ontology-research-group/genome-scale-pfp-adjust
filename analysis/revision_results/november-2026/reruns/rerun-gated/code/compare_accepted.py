"""Compare rerun outputs with accepted rev5 outputs (native-tau modes) and stage TSVs."""
import argparse, gzip, hashlib, json
from collections import Counter, defaultdict
from pathlib import Path

KEYS = ['n_present_complexes', 'n_singleton_complexes', 'complex_coherence', 'active_predictions',
        'promotions', 'demotions', 'process_coherence', 'known_taxon_consistent', 'gaef_taxon_consistent']
MODES = ['before', 'stage1_unknown', 'published_funtaxis', 'stage1', 'stage2_standard', 'stage2_demotion']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


ap = argparse.ArgumentParser()
ap.add_argument('--rerun', required=True); ap.add_argument('--accepted', required=True)
ap.add_argument('--taxa', required=True); ap.add_argument('--out', required=True)
a = ap.parse_args()
taxa = Path(a.taxa).read_text().split()
report = {}
for pred in ['seq_sim_timeset', 'deepgo-se', 'sprof-go']:
    c = Counter(); diffs = []
    for t in taxa:
        r = json.load(gzip.open(f'{a.rerun}/{pred}/cc/{t}.json.gz', 'rt'))
        x = json.load(gzip.open(f'{a.accepted}/{pred}/cc/{t}.json.gz', 'rt'))
        nat = x['native_tau']
        # input checksums by basename
        xi = {Path(k).name: v for k, v in x['inputs'].items()}
        for k, v in r['inputs'].items():
            if Path(k).name in xi and xi[Path(k).name] != v:
                diffs.append([t, 'input', Path(k).name])
        cx = x['code']; cr = r['code']
        c['code_files_differ'] = max(c['code_files_differ'], sum(cx.get(f) != h for f, h in cr.items()))
        rr = {(q['mode'], q['tau']): q for q in r['runs']}
        xr = {(q['mode'], q['tau']): q for q in x['runs']}
        for m in MODES:
            q, y = rr.get((m, nat)), xr.get((m, nat))
            if q is None or y is None:
                diffs.append([t, m, 'missing']); continue
            same = True
            for k in KEYS:
                if q.get(k) != y.get(k):
                    same = False; diffs.append([t, m, k, y.get(k), q.get(k)])
            if q['test_scores'] != y['test_scores']:
                same = False; diffs.append([t, m, 'test_scores'])
            if 'solver' in y:
                c[f'{m}_status_acc_{y["solver"]["status"]}'] += 1
                c[f'{m}_status_rerun_{q["solver"]["status"]}'] += 1
                if abs(y['solver']['objective'] - q['solver']['objective']) > 1e-6:
                    c[f'{m}_objective_differs'] += 1
            c[f'{m}_identical_metrics'] += same
        for m in ['stage1', 'stage2_standard', 'stage2_demotion']:
            fr = sorted(Path(f'{a.rerun}/{pred}/cc/{m}').glob(f'optimized_fold_*_taxon_{t}.tsv'))
            fx = sorted(Path(f'{a.accepted}/{pred}/cc/{m}').glob(f'optimized_fold_*_taxon_{t}.tsv'))
            ok = len(fr) == len(fx) == 1 and sha(fr[0]) == sha(fx[0])
            c[f'{m}_tsv_identical'] += ok
            if not ok:
                diffs.append([t, m, 'tsv_differs'])
    report[pred] = {'counts': dict(c), 'diffs': diffs}
    print(pred, dict(c), len(diffs))
Path(a.out).write_text(json.dumps(report, indent=1))
