"""How are Stage 1 demotions exported in the accepted outputs? (absent vs capped at tau)"""
import gzip, json, sys
from collections import Counter
from pathlib import Path
acc = sys.argv[1]
for pred in ['deepgo-se', 'seq_sim_timeset', 'sprof-go']:
    for asp in ['cc', 'mf', 'bp']:
        c = Counter()
        for f in sorted(Path(acc, pred, asp).glob('*.json.gz')):
            d = json.load(gzip.open(f, 'rt'))
            tau = d['native_tau']
            by = {r['mode']: r for r in d['runs'] if r['tau'] == tau}
            b, s = by['before']['test_scores'], by['stage1']['test_scores']
            for p, sc in b.items():
                o = s.get(p, {})
                for g, v in sc.items():
                    if g not in o:
                        c['absent_above' if v > tau else 'absent_below'] += 1
                    elif o[g] != v:
                        c['capped' if o[g] == min(v, tau) else 'other'] += 1
            c['active_before'] += by['stage1']['before']['active_predictions']
            c['demotions_all_proteins'] += by['stage1']['demotions']
            c['n_test'] += len(b)
            c['modified_test'] += len(by['stage1']['modified_test_proteins'])
        print(pred, asp, dict(c), flush=True)
