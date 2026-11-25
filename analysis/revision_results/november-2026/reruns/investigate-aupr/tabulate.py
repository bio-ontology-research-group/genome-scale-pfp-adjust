import json, sys, glob
from pathlib import Path
rows = []
for f in sorted(glob.glob(sys.argv[1] + '/filter_*.json')):
    d = json.load(open(f)); m = d['metrics']; s = d['stats']; ms = d['modified_subset']
    cc = d.get('crosscheck', {})
    ok = all(abs(cc[k]['aupr'] - m[k]['aupr']) < 1e-12 and abs(cc[k]['fmax'] - m[k]['fmax']) < 1e-12 for k in ['before', 'stage1'] if k in cc)
    def fmt(k):
        x = m.get(k)
        return '%.4f/%.3f/%.4f/%.4f' % (x['fmax'], x['smin'], x['avg_auc'], x['aupr']) if x else '-'
    ft = cc.get('published_funtaxis')
    print(f"{d['predictor']:16s} {d['aspect']} n={d['n_valid_proteins']:5d} check={ok} | before {fmt('before')} | stage1 {fmt('stage1')} (mod {ms['stage1']['n']}) | filter {fmt('filter_all')} (mod {ms['filter_all']['n']}) | unknown {fmt('stage1_unknown')} | dAUPR s1 {m['stage1']['aupr']-m['before']['aupr']:+.5f} filt {m['filter_all']['aupr']-m['before']['aupr']:+.5f} funtaxis {(ft['aupr']-m['before']['aupr']) if ft else float('nan'):+.5f} | removed {s.get('filter_removed',0)} (>tau {s.get('filter_removed_above_tau',0)}, TP {s.get('filter_removed_true_label',0)}) s1chg {s.get('stage1_changed',0)}")
    sub = ms['filter_all']
    if 'before' in sub and sub['before']:
        b, a = sub['before'], sub['after']
        print(f"    filter-modified subset n={sub['n']} taxa={sub['n_taxa']}: Fmax {b['fmax']:.4f}->{a['fmax']:.4f} Smin {b['smin']:.3f}->{a['smin']:.3f} AUC {b['avg_auc'] or 0:.4f}->{a['avg_auc'] or 0:.4f} AUPR {b['aupr']:.4f}->{a['aupr']:.4f}")
