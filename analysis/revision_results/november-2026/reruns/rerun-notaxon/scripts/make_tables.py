"""Tables for the no-taxon siglogit rerun: flips, consistency/GAEF, per-protein CAFA metrics, modified proteins."""
import json, sys
from collections import Counter, defaultdict
from pathlib import Path
new, acc, perorg, outdir = map(Path, sys.argv[1:5])
outdir.mkdir(parents=True, exist_ok=True)
NAMES = {'seq_sim_timeset': 'Seq-Sim', 'deepgo-se': 'DeepGO-SE', 'sprof-go': 'SPROF-GO', 'seq_sim_heldout': 'Seq-Sim', 'mlp_heldout': 'MLP-ESM2'}
SPLITS = {'timeset': ['seq_sim_timeset', 'deepgo-se', 'sprof-go'], 'heldout': ['seq_sim_heldout', 'mlp_heldout']}
MODES = ['before', 'stage1', 'stage1_unknown', 'stage1_unknown_siglogit', 'stage1_unknown_siglogit_aspect']
MLAB = {'before': 'Preds', 'stage1': 'Taxons', 'stage1_unknown': 'No taxons (current)',
        'stage1_unknown_siglogit': 'No taxons (siglogit, joint)', 'stage1_unknown_siglogit_aspect': 'No taxons (siglogit, per aspect)'}
def rows(d, pred, asp):
    f = d / f'{pred}_{asp}.json'
    if not f.exists():
        return None
    out = {}
    for r in json.loads(f.read_text())['rows']:
        tau_native = .3 if pred == 'deepgo-se' else .001
        if abs(r['tau'] - tau_native) < 1e-9:
            out[r['mode']] = r
    return out
def n(x): return f'{x:,}'.replace(',', '{,}')
def p(x, k=3): return '--' if x is None else f'{x:.{k}f}'
res = {'reproduction': [], 'flips': [], 'consistency': [], 'metrics': []}
tsv = defaultdict(list)
tex = {}
lineage = json.loads(perorg.read_text()) if perorg.exists() else []
porg = defaultdict(list)
for r in lineage:
    porg[r['predictor'], r['aspect']].append(r)
for split, preds in SPLITS.items():
    fl, co, me = [], [], []
    for pred in preds:
        for asp in ['cc', 'mf', 'bp']:
            N, A = rows(new, pred, asp), rows(acc, pred, asp)
            if N is None:
                continue
            # reproduction vs accepted
            for mode in ['before', 'stage1', 'stage1_unknown']:
                a, b = A[mode], N[mode]
                keys = ['demotions', 'active_input', 'n_modified_test_proteins', 'known_taxon_consistency', 'gaef_taxon_consistency']
                mk = ['fmax', 'smin', 'aupr', 'avg_auc']
                diff = {k: (a[k], b[k]) for k in keys if a[k] != b[k]}
                diff.update({k: (a['metrics'][k], b['metrics'][k]) for k in mk if abs(a['metrics'][k] - b['metrics'][k]) > 1e-12})
                res['reproduction'].append({'predictor': pred, 'aspect': asp, 'mode': mode, 'identical': not diff, 'diff': diff})
            act = N['before']['active_input']
            label = NAMES[pred]
            frow = [label, asp.upper(), act] + sum(([N[m]['demotions'], N[m]['demotion_active_pct']] for m in MODES[1:]), [])
            fl.append(frow); res['flips'].append(dict(split=split, predictor=pred, aspect=asp, active=act,
                **{m: {'demotions': N[m]['demotions'], 'pct': N[m]['demotion_active_pct']} for m in MODES[1:]}))
            # per-organism lineage + consistency counts
            orgs = porg.get((pred, asp), [])
            cats = Counter(o['runs']['stage1_unknown_siglogit']['lineage_inference']['category'] for o in orgs)
            catsa = Counter(o['runs']['stage1_unknown_siglogit_aspect']['lineage_inference']['category'] for o in orgs)
            for m in MODES:
                r = N[m]
                co.append([label, asp.upper(), MLAB[m], r['known_taxon_consistency'], r['gaef_taxon_consistency'], r['process_coherence'],
                           r['pathway_coherence'], r['essential_completeness'], r['complex_coherence']])
                me.append([label, asp.upper(), MLAB[m], r['metrics']['fmax'], r['metrics']['smin'], r['metrics']['aupr'], r['metrics']['avg_auc'],
                           r['n_modified_test_proteins'], r['n_annotated_proteins']])
                res['consistency'].append(dict(split=split, predictor=pred, aspect=asp, mode=m, n_organisms=r['n_organisms'],
                    known_taxon_consistency=r['known_taxon_consistency'], gaef_taxon_consistency=r['gaef_taxon_consistency'],
                    process_coherence=r['process_coherence'], pathway_coherence=r['pathway_coherence'],
                    essential_completeness=r['essential_completeness'], complex_coherence=r['complex_coherence'],
                    lineage_categories_joint=dict(cats), lineage_categories_aspect=dict(catsa)))
                res['metrics'].append(dict(split=split, predictor=pred, aspect=asp, mode=m, **{k: r['metrics'][k] for k in ['fmax', 'smin', 'aupr', 'avg_auc', 'precision', 'recall']},
                    n_modified_test_proteins=r['n_modified_test_proteins'], n_annotated_proteins=r['n_annotated_proteins'],
                    modified_metrics=r.get('modified_metrics') and {k: r['modified_metrics'][k] for k in ['fmax', 'smin', 'aupr', 'avg_auc']}))
    if not fl:
        continue
    # Flip table (Table 1 layout + new columns)
    L = [r'\begin{tabular}{llrrrrrrrrr}', r'\toprule',
         r'Model & Sub & Preds & \multicolumn{2}{c}{With taxons} & \multicolumn{2}{c}{No taxons (current)} & \multicolumn{2}{c}{No taxons (siglogit)} & \multicolumn{2}{c}{Siglogit per aspect} \\',
         r'\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}\cmidrule(lr){10-11}', r'& & & Flips & \% & Flips & \% & Flips & \% & Flips & \% \\', r'\midrule']
    last = None
    for r in fl:
        if last is not None and r[0] != last: L.append(r'\midrule')
        L.append(' & '.join([r[0] if r[0] != last else '', r[1], n(r[2])] + [n(r[i]) if i % 2 == 1 else p(r[i]) for i in range(3, 11)]) + r' \\')
        last = r[0]
    L += [r'\bottomrule', r'\end{tabular}']
    (outdir / f'flips_{split}.tex').write_text('% Stage 1 demotions of active pairs at native tau, summed over organisms; % = demotions / active input pairs\n' + '\n'.join(L) + '\n')
    def simple(rowsx, header, fmt, name, comment):
        L = [r'\begin{tabular}{lll' + 'r' * (len(header) - 3) + '}', r'\toprule', ' & '.join(header) + r' \\', r'\midrule']
        last = None
        for r in rowsx:
            key = (r[0], r[1])
            if last is not None and key != last: L.append(r'\midrule' if r[0] != last[0] else r'\addlinespace')
            L.append(' & '.join([r[0] if key != last else '', r[1] if key != last else '', r[2]] + [f(v) for f, v in zip(fmt, r[3:])]) + r' \\')
            last = key
        L += [r'\bottomrule', r'\end{tabular}']
        (outdir / f'{name}_{split}.tex').write_text(f'% {comment}\n' + '\n'.join(L) + '\n')
    pc = lambda v: '--' if v is None else f'{100 * v:.1f}'
    simple(co, ['Model', 'Sub', 'Setting', r'KTC \%', r'GTC \%', r'Proc.\ coh.\ \%', r'Path.\ coh.\ \%', r'Ess.\ compl.\ \%', r'Cplx.\ coh.\ \%'],
           [pc] * 6, 'consistency', 'KTC: repaired proteome consistent with the TRUE lineage (% organisms); GTC: GAEF pairwise screen; GAEF metrics macro-averaged over organisms')
    simple(me, ['Model', 'Sub', 'Setting', r'$F_{\max}$', r'$S_{\min}$', 'AUPR', 'AUC', 'Changed', 'Test'],
           [lambda v: f'{v:.4f}', lambda v: f'{v:.3f}', lambda v: f'{v:.4f}', lambda v: f'{v:.4f}', lambda v: n(v), lambda v: n(v)],
           'metrics', 'Per-protein CAFA metrics (aggregate_revision_audit evaluator); Changed = test proteins whose stored scores changed')
# Per-organism lineage breakdown: inference category vs true-lineage consistency after repair
res['lineage'] = []
for split, preds in SPLITS.items():
    L = [r'\begin{tabular}{llrrrrrrrr}', r'\toprule',
         r'Model & Sub & \multicolumn{3}{c}{Joint inference (C/W/U)} & \multicolumn{3}{c}{KTC (organisms)} & \multicolumn{2}{c}{Siglogit KTC by inference} \\',
         r'\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-10}',
         r'& & C & W & U & Current & Siglogit & Per aspect & C & W/U \\', r'\midrule']
    any_rows = False; last = None
    for pred in preds:
        for asp in ['cc', 'mf', 'bp']:
            orgs = porg.get((pred, asp))
            if not orgs:
                continue
            any_rows = True
            cat = lambda o: o['runs']['stage1_unknown_siglogit']['lineage_inference']['category']
            k = lambda o, m: bool(o['runs'][m]['known_taxon_consistent'])
            c = Counter(cat(o) for o in orgs)
            ktc = {m: sum(k(o, m) for o in orgs) for m in MODES[2:]}
            kc = sum(k(o, 'stage1_unknown_siglogit') for o in orgs if cat(o) == 'correct')
            kw = sum(k(o, 'stage1_unknown_siglogit') for o in orgs if cat(o) != 'correct')
            same = sum(o['runs']['stage1_unknown_siglogit']['demotions'] == o['runs']['stage1']['demotions'] for o in orgs)
            res['lineage'].append(dict(split=split, predictor=pred, aspect=asp, n=len(orgs), categories=dict(c), ktc=ktc,
                                       ktc_siglogit_correct=kc, ktc_siglogit_not_correct=kw, siglogit_equals_known_demotions=same,
                                       fixed_rank_counts=dict(Counter(o['runs']['stage1_unknown_siglogit']['lineage_inference']['deepest'] for o in orgs))))
            if last is not None and pred != last: L.append(r'\midrule')
            L.append(' & '.join([NAMES[pred] if pred != last else '', asp.upper(), str(c['correct']), str(c['wrong']), str(c['unresolved']),
                                 str(ktc['stage1_unknown']), str(ktc['stage1_unknown_siglogit']), str(ktc['stage1_unknown_siglogit_aspect']),
                                 f"{kc}/{c['correct']}", f"{kw}/{len(orgs) - c['correct']}"]) + r' \\')
            last = pred
    L += [r'\bottomrule', r'\end{tabular}']
    if any_rows:
        (outdir / f'lineage_{split}.tex').write_text('% C/W/U: joint siglogit inference correct prefix / wrong / unresolved; KTC: organisms whose repaired proteome is consistent with the true lineage\n' + '\n'.join(L) + '\n')
(outdir / 'results.json').write_text(json.dumps(res, indent=1))
for r in res['lineage']:
    print('LIN', r['predictor'], r['aspect'], r['categories'], r['ktc'], 'C-ktc', r['ktc_siglogit_correct'], 'notC-ktc', r['ktc_siglogit_not_correct'], 'same-as-known', r['siglogit_equals_known_demotions'], r['fixed_rank_counts'])
for r in res['reproduction']:
    print('REPRO', r['predictor'], r['aspect'], r['mode'], 'identical' if r['identical'] else r['diff'])
