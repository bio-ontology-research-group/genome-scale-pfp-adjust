"""Tables (_full suffix) comparing no-taxon repairs: current / siglogit required-only / siglogit full selected lineage / known taxon.

Usage: make_tables_full.py SUMMARIES_OLD SUMMARIES_FULL PERORG_OLD PERORG_FULL OUTDIR
"""
import json, sys
from collections import Counter, defaultdict
from pathlib import Path
old, full, perorg_old, perorg_full, outdir = map(Path, sys.argv[1:6])
outdir.mkdir(parents=True, exist_ok=True)
NAMES = {'seq_sim_timeset': 'Seq-Sim', 'deepgo-se': 'DeepGO-SE', 'sprof-go': 'SPROF-GO', 'seq_sim_heldout': 'Seq-Sim', 'mlp_heldout': 'MLP-ESM2'}
SPLITS = {'timeset': ['seq_sim_timeset', 'deepgo-se', 'sprof-go'], 'heldout': ['seq_sim_heldout', 'mlp_heldout']}
FULL = 'stage1_unknown_siglogit_full'
MODES = ['before', 'stage1', 'stage1_unknown', 'stage1_unknown_siglogit', FULL]
MLAB = {'before': 'Preds', 'stage1': 'Taxons', 'stage1_unknown': 'No taxons (current)',
        'stage1_unknown_siglogit': 'No taxons (siglogit, required)', FULL: 'No taxons (siglogit, full lineage)'}
MK = ['fmax', 'smin', 'aupr', 'avg_auc']


def rows(d, pred, asp):
    f = d / f'{pred}_{asp}.json'
    if not f.exists():
        return None
    tau_native = .3 if pred == 'deepgo-se' else .001
    return {r['mode']: r for r in json.loads(f.read_text())['rows'] if abs(r['tau'] - tau_native) < 1e-9}


def n(x): return f'{x:,}'.replace(',', '{,}')
def p(x, k=3): return '--' if x is None else f'{x:.{k}f}'
pc = lambda v: '--' if v is None else f'{100 * v:.1f}'

porg = {}
for path in [perorg_old, perorg_full]:
    for r in json.loads(path.read_text()):
        porg.setdefault((r['predictor'], r['aspect'], r['taxon']), {}).update(r['runs'])

res = {'before_check': [], 'flips': [], 'consistency': [], 'metrics': [], 'lineage': [], 'witness': []}


def table(rowsx, header, fmt, name, split, comment, keycols=2):
    L = [r'\begin{tabular}{' + 'l' * (keycols + 1) + 'r' * (len(header) - keycols - 1) + '}', r'\toprule', ' & '.join(header) + r' \\', r'\midrule']
    last = None
    for r in rowsx:
        key = tuple(r[:keycols])
        if last is not None and key != last:
            L.append(r'\midrule' if r[0] != last[0] else r'\addlinespace')
        L.append(' & '.join([r[0] if key != last else '', r[1] if key != last else ''] + list(r[2:keycols + 1]) +
                            [f(v) for f, v in zip(fmt, r[keycols + 1:])]) + r' \\')
        last = key
    L += [r'\bottomrule', r'\end{tabular}']
    (outdir / f'{name}_{split}_full.tex').write_text(f'% {comment}\n' + '\n'.join(L) + '\n')


for split, preds in SPLITS.items():
    fl, co, me, li, wi = [], [], [], [], []
    for pred in preds:
        for asp in ['cc', 'mf', 'bp']:
            O, F = rows(old, pred, asp), rows(full, pred, asp)
            if F is None:
                continue
            # the full rerun recomputes 'before': must equal the old baseline exactly
            diff = {k: (O['before']['metrics'][k], F['before']['metrics'][k]) for k in MK if O['before']['metrics'][k] != F['before']['metrics'][k]}
            res['before_check'].append(dict(predictor=pred, aspect=asp, identical=not diff, diff=diff))
            N = {**O, FULL: F[FULL]}
            act = N['before']['active_input']
            label = NAMES[pred]
            fl.append([label, asp.upper(), act] + sum(([N[m]['demotions'], N[m]['demotion_active_pct']] for m in MODES[1:]), []))
            res['flips'].append(dict(split=split, predictor=pred, aspect=asp, active=act,
                                     **{m: {'demotions': N[m]['demotions'], 'pct': N[m]['demotion_active_pct']} for m in MODES[1:]}))
            for m in MODES:
                r = N[m]
                co.append([label, asp.upper(), MLAB[m], r['known_taxon_consistency'], r['gaef_taxon_consistency'], r['process_coherence'],
                           r['pathway_coherence'], r['essential_completeness'], r['complex_coherence']])
                me.append([label, asp.upper(), MLAB[m], *(r['metrics'][k] for k in MK), r['n_modified_test_proteins'], r['n_annotated_proteins']])
                res['consistency'].append(dict(split=split, predictor=pred, aspect=asp, mode=m, n_organisms=r['n_organisms'],
                                               **{k: r[k] for k in ['known_taxon_consistency', 'gaef_taxon_consistency', 'process_coherence',
                                                                    'pathway_coherence', 'essential_completeness', 'complex_coherence']}))
                res['metrics'].append(dict(split=split, predictor=pred, aspect=asp, mode=m, **{k: r['metrics'][k] for k in MK + ['precision', 'recall']},
                                           n_modified_test_proteins=r['n_modified_test_proteins'], n_annotated_proteins=r['n_annotated_proteins'],
                                           modified_metrics=r.get('modified_metrics') and {k: r['modified_metrics'][k] for k in MK}))
            orgs = [v for (pp, a, t), v in porg.items() if pp == pred and a == asp]
            orgs = [o for o in orgs if FULL in o]
            k = lambda o, m: bool(o[m]['known_taxon_consistent'])
            w = lambda o: o[FULL]['witness']
            ok = [o for o in orgs if w(o)['deepest_in_true_lineage']]
            bad = [o for o in orgs if w(o)['deepest_in_true_lineage'] is False]
            ktc = {m: sum(k(o, m) for o in orgs) for m in MODES[2:]}
            same = sum(o[FULL]['demotions'] == o['stage1']['demotions'] for o in orgs)
            li.append([label, asp.upper(), str(len(orgs)), ktc['stage1_unknown'], ktc['stage1_unknown_siglogit'], ktc[FULL],
                       f"{sum(k(o, FULL) for o in ok)}/{len(ok)}", f"{sum(k(o, FULL) for o in bad)}/{len(bad)}", same])
            res['lineage'].append(dict(split=split, predictor=pred, aspect=asp, n=len(orgs), ktc=ktc,
                                       ktc_full_deepest_correct=[sum(k(o, FULL) for o in ok), len(ok)],
                                       ktc_full_deepest_wrong=[sum(k(o, FULL) for o in bad), len(bad)],
                                       full_equals_known_demotions=same,
                                       full_modes=dict(Counter(o[FULL]['full_mode'] for o in orgs))))
            if asp == 'cc':
                W = [w(o) for o in orgs]
                c = dict(n=len(W), fallback=sum(x['fallback'] for x in W),
                         beyond=sum(bool(x['deepest_beyond_required']) for x in W),
                         beyond_ranked=sum(bool(x['deepest_ranked_beyond_required']) for x in W),
                         deepest_wrong=sum(x['deepest_in_true_lineage'] is False for x in W),
                         beyond_and_wrong=sum(bool(x['deepest_beyond_required']) and x['deepest_in_true_lineage'] is False for x in W),
                         beyond_and_correct=sum(bool(x['deepest_beyond_required']) and bool(x['deepest_in_true_lineage']) for x in W),
                         required_wrong=sum(x['saved_category'] == 'wrong' for x in W),
                         reproduces_saved=sum(bool(x['reproduces_saved']) for x in W),
                         deepest_ranked_rank=dict(Counter(x['deepest_ranked_rank'] for x in W)))
                res['witness'].append(dict(split=split, predictor=pred, **c))
                wi.append([label, c['n'], c['beyond'], c['beyond_ranked'], c['beyond_and_correct'], c['beyond_and_wrong'], c['deepest_wrong'], c['required_wrong'], c['fallback']])
    if not fl:
        continue
    L = [r'\begin{tabular}{llrrrrrrrrr}', r'\toprule',
         r'Model & Sub & Preds & \multicolumn{2}{c}{With taxons} & \multicolumn{2}{c}{No taxons (current)} & \multicolumn{2}{c}{Siglogit required} & \multicolumn{2}{c}{Siglogit full lineage} \\',
         r'\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}\cmidrule(lr){10-11}', r'& & & Flips & \% & Flips & \% & Flips & \% & Flips & \% \\', r'\midrule']
    last = None
    for r in fl:
        if last is not None and r[0] != last: L.append(r'\midrule')
        L.append(' & '.join([r[0] if r[0] != last else '', r[1], n(r[2])] + [n(r[i]) if i % 2 == 1 else p(r[i]) for i in range(3, 11)]) + r' \\')
        last = r[0]
    L += [r'\bottomrule', r'\end{tabular}']
    (outdir / f'flips_{split}_full.tex').write_text('% Stage 1 demotions of active pairs at native tau, summed over organisms; % = demotions / active input pairs\n' + '\n'.join(L) + '\n')
    table(co, ['Model', 'Sub', 'Setting', r'KTC \%', r'GTC \%', r'Proc.\ coh.\ \%', r'Path.\ coh.\ \%', r'Ess.\ compl.\ \%', r'Cplx.\ coh.\ \%'],
          [pc] * 6, 'consistency', split, 'KTC: repaired proteome consistent with the TRUE lineage (% organisms); GTC: GAEF pairwise screen; GAEF metrics macro-averaged over organisms')
    table(me, ['Model', 'Sub', 'Setting', r'$F_{\max}$', r'$S_{\min}$', 'AUPR', 'AUC', 'Changed', 'Test'],
          [lambda v: f'{v:.4f}', lambda v: f'{v:.3f}', lambda v: f'{v:.4f}', lambda v: f'{v:.4f}', n, n],
          'metrics', split, 'Per-protein CAFA metrics (aggregate_revision_audit evaluator); Changed = test proteins whose stored scores changed')
    table(li, ['Model', 'Sub', 'Organisms', 'Current', 'Required', 'Full', 'Full | sel.\\ correct', 'Full | sel.\\ wrong', 'Full = known'],
          [str, str, str, str, str, str, str], 'lineage', split,
          'Organisms whose repaired proteome is consistent with the true lineage (KTC counts); sel. correct/wrong: deepest taxon of the selected siglogit lineage in / not in the true lineage; Full = known: identical demotion count to known-taxon Stage 1')
    L = [r'\begin{tabular}{lrrrrrrrr}', r'\toprule',
         r'Model & Organisms & Beyond req. & Beyond req.\ (ranked) & Beyond, correct & Beyond, wrong & Selected wrong & Required wrong & Fallback \\', r'\midrule']
    L += [' & '.join([r[0]] + [str(v) for v in r[1:]]) + r' \\' for r in wi]
    L += [r'\bottomrule', r'\end{tabular}']
    (outdir / f'witness_{split}_full.tex').write_text(
        '% Joint siglogit selected optimum (first CP-SAT solve: seed 42, 1 worker, sorted variables). Beyond req.: deepest selected taxon is not certified (required);'
        ' ranked: deepest selected RANKS_7 taxon deeper than deepest required one; wrong: deepest selected taxon not in the true lineage;'
        ' Required wrong: certified prefix itself wrong; Fallback: no ranked taxon selected (current no-taxon repair)\n' + '\n'.join(L) + '\n')
(outdir / 'results_full.json').write_text(json.dumps(res, indent=1))
for r in res['before_check']:
    print('BEFORE', r['predictor'], r['aspect'], 'identical' if r['identical'] else r['diff'])
for r in res['witness']:
    print('WIT', r)
for r in res['lineage']:
    print('LIN', r['predictor'], r['aspect'], r['ktc'], r['ktc_full_deepest_correct'], r['ktc_full_deepest_wrong'], 'same-as-known', r['full_equals_known_demotions'], r['full_modes'])
for r in res['metrics']:
    print('MET', r['split'], r['predictor'], r['aspect'], r['mode'], *(f"{r[k]:.4f}" for k in MK))
