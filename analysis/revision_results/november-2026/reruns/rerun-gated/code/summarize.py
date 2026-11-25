"""Summarize the gated rerun into JSON + generated-table-style LaTeX snippets (not for direct inclusion)."""
import gzip, json, re, statistics
from collections import Counter, defaultdict
from pathlib import Path

R = Path('/mnt/data1/leechuck/genoadjust-ws/rerun-gated')
PREDS = [('seq_sim_timeset', 'Seq-Sim', 'SS-E', .001), ('deepgo-se', 'DeepGO-SE', 'DGO-E', .3), ('sprof-go', 'SPROF-GO', 'SP-E', .001)]
MODES = ['before', 'stage1', 'stage2_standard', 'stage2_demotion', 'stage2_gated']
taxa = (R/'taxa.txt').read_text().split()
res = {}
for pred, name, short, tau in PREDS:
    P = res[pred] = {}
    for src, path in [('rerun', R/'summaries'/f'{pred}_cc.json'), ('accepted', R/'accepted-summaries'/f'{pred}_cc.json')]:
        rows = {r['mode']: r for r in json.loads(path.read_text())['rows'] if abs(r['tau']-tau) < 1e-12}
        P[src] = {m: {k: rows[m][k] for k in ['complex_coherence', 'process_coherence', 'known_taxon_consistency', 'gaef_taxon_consistency',
                                               'promotions', 'demotions', 'solver_status_counts', 'solver_objective_sum', 'solver_bound_sum',
                                               'n_modified_test_proteins', 'n_annotated_proteins']}
                  | {'metrics': {k: rows[m]['metrics'][k] for k in ['fmax', 'smin', 'aupr', 'avg_auc']},
                     'modified_metrics': rows[m]['modified_metrics'] and {k: rows[m]['modified_metrics'][k] for k in ['fmax', 'smin', 'aupr', 'avg_auc']}}
                  for m in MODES if m in rows}
    org = defaultdict(lambda: defaultdict(list))
    gated = Counter()
    for t in taxa:
        d = json.load(gzip.open(R/'out'/pred/'cc'/f'{t}.json.gz', 'rt'))
        for q in d['runs']:
            if q['mode'] in MODES[2:]:
                o = org[q['mode']]
                o['edits'].append(q['promotions'] + q['demotions'])
                o['seconds'].append(q['solver']['seconds'])
                o['status'].append(q['solver']['status'])
                o['singletons_after'].append(q['n_singleton_complexes'])
                o['taxon_ok'].append(q['known_taxon_consistent'])
                if q['mode'] == 'stage2_gated':
                    o['formulation'].append(q['solver'].get('formulation'))
                    gated['singletons'] += q['gated']['n_singletons']
                    gated['with_evidence'] += q['gated']['n_with_evidence']
                    gated.update('outcome_'+k for k, v in q['gated']['outcome_counts'].items() for _ in range(v))
            if q['mode'] == 'stage1':
                org['stage1']['seconds'].append(q['seconds_stage1'])
    wall = []
    for t in taxa:
        m = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', (R/'logs'/f'{pred}_{t}.log').read_text())
        parts = [float(x) for x in m.group(1).split(':')]
        wall.append(sum(v * 60**i for i, v in enumerate(reversed(parts))))
    P['per_organism'] = {m: {'organisms_with_edits': sum(e > 0 for e in o['edits']) if 'edits' in o else None,
                             'status': dict(Counter(o['status'])), 'formulation': dict(Counter(o['formulation'])),
                             'solver_seconds_sum': sum(o['seconds']), 'solver_seconds_median': statistics.median(o['seconds']),
                             'solver_seconds_max': max(o['seconds']),
                             'organisms_with_singletons_after': sum(s > 0 for s in o['singletons_after']),
                             'taxon_consistent_organisms': sum(o['taxon_ok'])} if m != 'stage1' else
                         {'seconds_sum': sum(o['seconds']), 'seconds_max': max(o['seconds'])} for m, o in org.items()}
    P['gated_counts'] = dict(gated)
    P['job_wall_seconds'] = {'sum': sum(wall), 'median': statistics.median(wall), 'max': max(wall)}
    P['membership'] = {}
    for m in ['standard', 'demotion', 'gated']:
        x = json.loads((R/'membership'/f'membership_{pred}_{m}.json').read_text())
        P['membership'][m] = {'promotions': x['promotions'], 'demotions': x['demotions'],
                              'unmapped_p': x['unmapped_promotions_with_reference'], 'unmapped_d': x['unmapped_demotions_with_reference']}
ev = json.loads((R/'eval'/'eval_gated.json').read_text())
for pred, *_ in PREDS:
    v = ev[pred]['sx_prefer']
    res[pred]['uniprot_cp'] = {k: v.get(k, 0) for k in ['initial_singletons', 'outcome_promoted', 'outcome_demoted', 'outcome_unresolved', 'X_prom', 'X_prom_uni',
                                                        'X_prom_cp_assess', 'X_prom_cp_listed', 'X_dem', 'X_dem_cp_assess', 'X_dem_cp_listed', 'X_dem_lone_uni',
                                                        'prom_cx', 'prom_cx', 'prom_uni_any', 'prom_cp_assess', 'prom_cp_listed', 'dem_cx', 'dem_cp_assess', 'dem_cp_listed',
                                                        'prom_test', 'prom_test_inref', 'dem_test', 'dem_test_inref', 'optimal']}
(R/'summaries'/'gated_summary.json').write_text(json.dumps(res, indent=1))

# ---- LaTeX snippets ----
def pct(x, bold=True):
    if x is None: return '--'
    v = 100*x
    return r'\textbf{100.0}' if bold and abs(v-100) < 1e-9 else f'{v:.2f}'
T = R/'tables'; T.mkdir(exist_ok=True)
lines = ['% Rerun (ws, gated driver) CC columns of tab:gaef_preds: Process% & Complex% & Taxon%.',
         '% adjusted$^\\dagger$ = standard Stage 2; adjusted$^\\ddagger$ = evidence-gated Stage 2 (sx_prefer). Not inserted into the paper.']
for pred, name, short, tau in PREDS:
    r = res[pred]['rerun']
    lines.append(f'% {name}')
    for m, lab in [('before', 'preds'), ('stage1', 'adjusted'), ('stage2_standard', r'adjusted$^\dagger$'), ('stage2_gated', r'adjusted$^\ddagger$')]:
        x = r[m]
        lines.append(f'   & {lab:<20} & {pct(x["process_coherence"], False)} & {pct(x["complex_coherence"]) if m.startswith("stage2") else pct(x["complex_coherence"], False)} & '
                     f'{pct(x["known_taxon_consistency"]) if m != "before" else pct(x["known_taxon_consistency"], False)} \\\\')
(T/'gaef-cc-rows.tex').write_text('\n'.join(lines)+'\n')

lab = {'stage1': 'Stage 1', 'stage2_standard': 'Stage 2 standard', 'stage2_demotion': 'Stage 2 demotion', 'stage2_gated': 'Stage 2 gated'}
out = [r'\begin{tabular}{llrrrrr}', r'\toprule', r'Model & Mode & $F_{\max}$ & $S_{\min}$ & AUPR & AUC & Modified \\', r'\midrule']
for pred, name, short, tau in PREDS:
    for m in ['before', 'stage1', 'stage2_standard', 'stage2_demotion', 'stage2_gated']:
        x = res[pred]['rerun'][m]; k = x['metrics']
        out.append(f'{short} & {lab.get(m, "Predictions")} & {k["fmax"]:.4f} & {k["smin"]:.4f} & {k["aupr"]:.4f} & {k["avg_auc"]:.4f} & {x["n_modified_test_proteins"]}/{x["n_annotated_proteins"]} \\\\')
out += [r'\bottomrule', r'\end{tabular}']
(T/'metrics-gated.tex').write_text('% CC per-protein metrics on the timeset test proteins (paper evaluator). Modified: test proteins changed vs the mode input.\n' + '\n'.join(out)+'\n')

out = [r'\begin{tabular}{lrrrrrr}', r'\toprule', r'Model & Mode & Optimal & Feasible & Objective sum & Bound sum & Edited taxa \\', r'\midrule']
for pred, name, short, tau in PREDS:
    for m in ['stage2_standard', 'stage2_demotion', 'stage2_gated']:
        x = res[pred]['rerun'][m]; s = x['solver_status_counts']
        out.append(f'{short} & {lab[m]} & {s.get("OPTIMAL", 0)} & {s.get("FEASIBLE", 0)} & {x["solver_objective_sum"]:.4f} & {x["solver_bound_sum"]:.4f} & {res[pred]["per_organism"][m]["organisms_with_edits"]} \\\\')
out += [r'\bottomrule', r'\end{tabular}']
(T/'supp-solvers-gated.tex').write_text('% Table S17-style rows incl. gated (gated: 3 CP-SAT workers, 120 s).\n' + '\n'.join(out)+'\n')

out = [r'\begin{tabular}{llrrrr}', r'\toprule', r'Model & Mode & Promotions & Demotions & Unmapped P & Unmapped D \\', r'\midrule']
for pred, name, short, tau in PREDS:
    for m in ['standard', 'demotion', 'gated']:
        x = res[pred]['membership'][m]; p, d = x['promotions'], x['demotions']
        out.append(f'{short} & {m} & {p["total"]}/{p["with_curated_reference"]}/{p["true_subunit_promotions"]} & {d["total"]}/{d["with_curated_reference"]}/{d["curated_subunit_demotions"]} & {x["unmapped_p"]} & {x["unmapped_d"]} \\\\')
out += [r'\bottomrule', r'\end{tabular}']
(T/'supp-membership-gated.tex').write_text('% Table S18-style rows (total / assessable / listed subunit).\n' + '\n'.join(out)+'\n')
print(json.dumps({p: {'table': {m: [res[p]['rerun'][m][k] for k in ['process_coherence', 'complex_coherence', 'known_taxon_consistency']] for m in MODES},
                      'per_org': res[p]['per_organism'], 'gated': res[p]['gated_counts'], 'uni': res[p]['uniprot_cp'], 'wall': res[p]['job_wall_seconds']} for p in res}, indent=1))
