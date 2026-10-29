#!/usr/bin/env python3
"""Generate revision tables from a complete, validated result bundle.

No manuscript tables or rendered figures are used as numerical inputs.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

MODELS = ['seq_sim_timeset', 'deepgo-se', 'sprof-go', 'seq_sim_heldout', 'mlp_heldout']
SHORT = dict(zip(MODELS, ['SS-E', 'DGO-E', 'SP-E', 'SS-H', 'MLP-H']))
ASPECTS = ['cc', 'mf', 'bp']
RANKS = ['none', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
LEGEND = 'SS: Seq-Sim; DGO: DeepGO-SE; SP: SPROF-GO; MLP: MLP-ESM2. E: entry-date split; H: heldout split.'
MODES = {'before': 'Raw', 'stage1': 'Known taxon', 'stage1_unknown': 'Inferred taxon',
         'published_funtaxis': 'FunTaxIS-lite', 'stage2_standard': 'Stage 2 standard',
         'stage2_demotion': 'Stage 2 demotion'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fmt(value, digits=4):
    if value is None:
        return '--'
    return f'{value:.{digits}f}'


def pct(value):
    return fmt(100*value, 2) if value is not None else '--'


def compact(value):
    return '--' if value is None else f'{value:.4g}'


def mode_name(mode):
    if mode.startswith('soft_stage'):
        stage, penalty = mode.removeprefix('soft_stage').split('_')
        return f'Soft {stage}, {penalty}'
    return MODES[mode]


def render_table(out, name, caption, headers, rows, main=False):
    # The TSV and TeX are generated from the same in-memory rows.
    with (out/(name+'.tsv')).open('w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(rows)
    columns = 'l' + 'r'*(len(headers)-1)
    if main:
        lines = [r'\begin{table*}[t]', r'\centering\small', r'\caption{'+caption+'}',
                 r'\label{tab:'+name+'}', r'\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{'+columns+'}']
    else:
        lines = [r'\begingroup\scriptsize', r'\setlength{\tabcolsep}{4pt}', r'\begin{longtable}{'+columns+'}',
                 r'\caption{'+caption+r'}\label{tab:'+name+r'}\\']
    heading = ' & '.join(headers)+r' \\'
    lines += [r'\toprule', heading, r'\midrule']
    if not main:
        lines += [r'\endfirsthead', r'\toprule', heading, r'\midrule', r'\endhead']
    lines += [' & '.join(map(str,row))+r' \\' for row in rows]
    lines += [r'\bottomrule', r'\end{tabular}' if main else r'\end{longtable}', r'\end{table*}' if main else r'\endgroup']
    (out/(name+'.tex')).write_text('\n'.join(lines)+'\n')


def generate(args):
    source, out = Path(args.results), Path(args.output)
    required = [source/'summaries'/f'{m}_{a}.json' for m in MODELS for a in ASPECTS]
    required += [source/'lineage.json', source/'evidence-sensitivity.json', source/'export-validation.json']
    required += [source/f'membership_{m}_{mode}.json' for m in MODELS[:3] for mode in ['standard','demotion']]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError(f'Incomplete result bundle ({len(missing)} missing): {missing[:5]}')
    rows, index, sources = [], {}, {}
    code = None
    ic_inputs = None
    for path in required[:15]:
        data = json.loads(path.read_text())
        if data.get('schema_version') != 2 or not data['sources']:
            raise ValueError(f'Unaccepted summary schema: {path}')
        if not data.get('reference_manifest_sha256'):
            raise ValueError('Accuracy summaries require source-rebuilt reference labels')
        current_ic = (data['ic_reference_inputs'], data['ic_split_inputs'], data['reference_manifest_sha256'])
        if ic_inputs is None:
            ic_inputs = current_ic
        elif current_ic != ic_inputs:
            raise ValueError('Summaries use different IC reference inputs')
        for item in data['sources']:
            sources[tuple(item['task'])] = item
            if code is None:
                code = item['code']
            elif item['code'] != code:
                raise ValueError('Mixed experiment code manifests')
        for row in data['rows']:
            key = row['predictor'], row['aspect'], row['mode'], row['tau']
            if key in index:
                raise ValueError(f'Duplicate summary condition {key}')
            if row['n_organisms'] != (219 if key[0].endswith('_heldout') else 50):
                raise ValueError(f'Unexpected organism count: {key}')
            index[key] = row
            rows.append(row)
    def get(m, a, mode, tau=None):
        return index[m, a, mode, (.3 if m == 'deepgo-se' else .001) if tau is None else tau]
    for m in MODELS:
        for a in ASPECTS:
            for mode in ['before','stage1','stage1_unknown']:
                get(m,a,mode)
            if m != 'seq_sim_heldout':
                for t in [.001,.01,.05,.1,.2,.3,.4,.5]:
                    get(m,a,'before',t)
                    get(m,a,'stage1',t)
    export_qc = json.loads((source/'export-validation.json').read_text())
    if export_qc.get('validated') != 300 or len(export_qc.get('records', [])) != 300:
        raise ValueError('Missing independent hard-repair export validation')
    lineage = json.loads((source/'lineage.json').read_text())
    evidence = json.loads((source/'evidence-sensitivity.json').read_text())
    if len(evidence['rows']) != 50 or len(lineage['summary']) != 20:
        raise ValueError('Incomplete evidence or lineage summary')
    if lineage['input_sha256'] != {'/'.join(task): item['result_sha256'] for task, item in sources.items()}:
        raise ValueError('Lineage and metric summaries use different result files')
    checks = {(r['predictor'], r['taxon'], r['mode']): r for r in export_qc['records']}
    expected_checks = {(m, t, mode) for m, a, t in sources if m in MODELS[:3] and a == 'cc'
                       for mode in ['standard', 'demotion']}
    if set(checks) != expected_checks or len(checks) != 300:
        raise ValueError('Export validation organism manifest differs')
    for (m, t, mode), check in checks.items():
        if check['metadata_sha256'] != sources[m, 'cc', t]['result_sha256']:
            raise ValueError('Export validation uses different result metadata')
    memberships = {}
    for m in MODELS[:3]:
        for mode in ['standard','demotion']:
            d = json.loads((source/f'membership_{m}_{mode}.json').read_text())
            if d['n_organism_pairs'] != 50:
                raise ValueError('Incomplete curated membership comparison')
            for pair in d['provenance']['pairs']:
                check = checks[m, pair['taxon'], mode]
                if pair['stage1'] != check['stage1_sha256'] or pair['stage2'] != check['stage2_sha256']:
                    raise ValueError('Membership and export checks use different prediction files')
            if {p['taxon'] for p in d['provenance']['pairs']} != {t for mm, t, md in checks if mm == m and md == mode}:
                raise ValueError('Membership organism manifest differs')
            for key in ['promotions', 'demotions']:
                total = sum(c[key] for (mm, _, md), c in checks.items() if mm == m and md == mode)
                if total != d[key]['total'] or total != get(m, 'cc', 'stage2_'+mode)[key]:
                    raise ValueError('Serialized, membership and summary edit counts disagree')
            memberships[m,mode] = d
    destination = out
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = tempfile.TemporaryDirectory(prefix='.revision-tables-', dir=destination.parent)
    out = Path(staging.name)
    flip_rows, metric_rows, scope_rows, accuracy_rows, identifier_rows = [], [], [], [], []
    for m in MODELS:
        for a in ASPECTS:
            r, b = get(m,a,'stage1'), get(m,a,'before')
            label = SHORT[m]+' '+a.upper()
            flip_rows.append([label, fmt(r['tau'],3), f"{r['active_input']:,}", f"{r['demotions']:,}", fmt(r['demotion_active_pct'],3), compact(r['flip_grid_pct'])])
            paired = lambda key: fmt(b['metrics'][key])+'/'+fmt(r['metrics'][key])
            metric_rows.append([label,r['n_annotated_proteins'],pct(b['known_taxon_consistency'])+'/'+pct(r['known_taxon_consistency']),paired('fmax'),paired('smin'),paired('aupr')])
            scope_rows.append([label,r['n_organisms'],r['scope_counts']['n_mask_proteins'],r['scope_counts']['n_annotation_matches'],r['n_annotated_proteins'],r['n_annotated_organisms']])
            identifier_rows.append([label]+[r['identifier_normalization'].get(key,0) for key in
                ['alternative_id_pairs','merged_duplicate_pairs','unknown_or_obsolete_pairs','other_namespace_pairs','hierarchy_added_pairs','hierarchy_raised_pairs']])
            for mode in ['before','stage1','stage1_unknown']:
                rr=get(m,a,mode)
                accuracy_rows.append([label,mode_name(mode)]+[fmt(rr['metrics'][k]) for k in ['fmax','smin','avg_auc','aupr','precision','recall']])
    render_table(out,'main-flips', 'Stage 1 edit burden with the known taxon at native thresholds. Active and grid rates use different, explicitly defined denominators. '+LEGEND,
        ['Model/aspect',r'$\tau$','Active input','Demotions','Active '+r'\%','Grid '+r'\%'],flip_rows,True)
    render_table(out,'main-metrics','Raw/known-taxon Stage 1 results at native thresholds. TC is known-taxon consistency across all designated organisms; $n$ is the eligible protein count for accuracy. '+LEGEND,
        ['Model/aspect','$n$','TC '+r'\%',r'$F_{\max}$',r'$S_{\min}$','AUPR'],metric_rows,True)
    render_table(out,'supp-scope','Evaluation populations. Matched entries can lack non-root labels in a particular aspect. '+LEGEND,
        ['Model/aspect','Organisms','Candidates','Matched','Eligible proteins','Eligible organisms'],scope_rows)
    render_table(out,'supp-identifiers','Input score alignment before every comparison. Alternative IDs are mapped to primary IDs, duplicate scores merged by maximum, and unknown/obsolete or out-of-aspect scores excluded. Ancestor scores are then raised to the largest descendant score. Added/raised pairs are preprocessing operations, not repair flips. Counts precede thresholding and are counted once per input file. '+LEGEND,
        ['Model/aspect','Aliases','Merged','Unknown/obsolete','Other aspect','Added ancestors','Raised scores'],identifier_rows)
    render_table(out,'supp-accuracy','Accuracy at native thresholds. AUC excludes terms without both label classes. '+LEGEND,
        ['Model/aspect','Mode',r'$F_{\max}$',r'$S_{\min}$','AUC','AUPR','Precision','Recall'],accuracy_rows)
    tflips,tmetrics,tstructure=[],[],[]
    for m in [MODELS[0],MODELS[1],MODELS[2],MODELS[4]]:
        for a in ASPECTS:
            for t in [.001,.01,.05,.1,.2,.3,.4,.5]:
                r=get(m,a,'stage1',t)
                label=SHORT[m]+' '+a.upper()
                tflips.append([label,fmt(t,3),f"{r['active_input']:,}",f"{r['demotions']:,}",fmt(r['demotion_active_pct'],3),compact(r['flip_grid_pct']),pct(r['known_taxon_consistency']),pct(r['gaef_taxon_consistency'])])
                tmetrics.append([label,fmt(t,3)]+[fmt(r['metrics'][k]) for k in ['fmax','smin','avg_auc','aupr']])
                b=get(m,a,'before',t)
                structural = []
                for key in ['complex_coherence','process_coherence','pathway_coherence','essential_completeness']:
                    values=[]
                    for item in [b,r]:
                        value=item[key]
                        n=item['n_organisms'] if key=='complex_coherence' and value is not None else item.get(key+'_n',0)
                        values.append(pct(value)+' ('+str(n)+')' if value is not None else '--')
                    structural.append(' / '.join(values))
                tstructure.append([label,fmt(t,3)]+structural)
    render_table(out,'supp-threshold-flips','Threshold sweep: edit counts and structural consistency after Stage 1. KTC uses the known taxon; GTC is the bundled GAEF pairwise screen. '+LEGEND,
        ['Model/aspect',r'$\tau$','Active input','Demotions','Active '+r'\%','Grid '+r'\%','KTC '+r'\%','GTC '+r'\%'],tflips)
    render_table(out,'supp-threshold-metrics','Threshold sweep: accuracy after Stage 1. Raw continuous-score baselines are in Table S2 and do not depend on the repair threshold. '+LEGEND,
        ['Model/aspect',r'$\tau$',r'$F_{\max}$',r'$S_{\min}$','AUC','AUPR'],tmetrics)
    render_table(out,'supp-threshold-structure','Other structural measures before / after known-taxon Stage 1, as percentages (assessable organisms). Complex coherence assigns 100 percent to an organism without a present complex; process and pathway means exclude organisms with no assessable relations or pathways. Essential completeness applies to BP. '+LEGEND,
        ['Model/aspect',r'$\tau$','Complex','Process','Pathway','Essential'],tstructure)
    lrows=[]
    for r in lineage['summary']:
        expected=219 if r['predictor'].endswith('_heldout') else 50
        if r['n'] != expected:
            raise ValueError('Incomplete lineage population')
        for measure,key in [('Deepest','deepest_correct_rank'),('Divergence','first_divergence_rank')]:
            counts=r['supported'][key]
            if sum(counts.values()) != expected:
                raise ValueError('Lineage histogram does not sum to population')
            lrows.append([SHORT[r['predictor']]+' '+r['aspect'].upper(),measure]+[counts.get(rank,0) for rank in RANKS])
    render_table(out,'supp-lineage','Ranks required by the selected retained GO set. Deepest: deepest correct rank; Divergence: first divergent rank. None denotes no resolved correct rank or no observed divergence, respectively, and is not itself evidence of correct identity. Raw assignments and cross-aspect rank records are included in lineage.json. '+LEGEND,
        ['Model/aspect','Measure','None','Domain','Phylum','Class','Order','Family','Genus','Species'],lrows)
    joint_rows = [r for r in lineage['summary'] if r['aspect'] == 'joint']
    cross = [r'\par Cross-aspect compatibility / all three aspects resolved identically, by rank.',
             r'\begin{center}\scriptsize\begin{tabular}{lrrrrrrr}\toprule',
             'Model & Domain & Phylum & Class & Order & Family & Genus & Species '+r'\\', r'\midrule']
    for row in joint_rows:
        values = [str(row['cross_aspect_rank_agreement'][rank]['compatible'])+'/'+str(row['cross_aspect_rank_agreement'][rank]['all_three_resolved_equal']) for rank in RANKS[1:]]
        cross.append(' & '.join([SHORT[row['predictor']]]+values)+r'\\')
    cross += [r'\bottomrule\end{tabular}\end{center}']
    with (out/'supp-lineage.tex').open('a') as fh:
        fh.write('\n'.join(cross)+'\n')
    comparisons=[]
    conditions=['before','stage1','published_funtaxis']+[f'soft_stage{s}_{p}' for p in [.2,.5,1.0] for s in [1,2]]+['stage2_standard','stage2_demotion']
    for m in MODELS[:3]:
        for aspect in ASPECTS:
            for mode in conditions if aspect == 'cc' else ['before','stage1','published_funtaxis']:
                r=get(m,aspect,mode)
                comparisons.append([SHORT[m]+' '+aspect.upper(),mode_name(mode),pct(r['known_taxon_consistency']),pct(r['complex_coherence'])]+[fmt(r['metrics'][k]) for k in ['fmax','smin','avg_auc','aupr']])
    render_table(out,'supp-comparators','Published filtering in all three GO aspects, with complex repair and soft-constraint ablation in CC, on all 50 entry-date organisms per predictor. KTC and complex coherence are percentages; accuracy uses the exact eligible proteins in Table S1. '+LEGEND,
        ['Model/aspect','Mode','KTC '+r'\%','Coherence '+r'\%',r'$F_{\max}$',r'$S_{\min}$','AUC','AUPR'],comparisons)
    erows=[]
    for r in evidence['rows']:
        ref=r['complex_counts']['curated']; all_=r['complex_counts']['all_proteome']
        erows.append([r['taxon'],ref['present'],ref['singletons'],r['curated_singletons_with_partner_all_same_proteins'],r['curated_singletons_with_partner_all_proteome'],all_['present'],all_['singletons']])
    render_table(out,'supp-evidence','Annotation evidence sensitivity. Ref: reference evidence including TAS/IC. Resolved-same and resolved-all count original reference singletons that gain another annotated protein using all evidence on the same eligible proteins or the whole proteome. All-present/all-single count the expanded complex vocabulary, not the original denominator.',
        ['Taxon','Ref-present','Ref-single','Resolved-same','Resolved-all','All-present','All-single'],erows)
    reference_singletons=sum(r['complex_counts']['curated']['singletons'] for r in evidence['rows'])
    same=sum(r['curated_singletons_with_partner_all_same_proteins'] for r in evidence['rows'])
    whole=sum(r['curated_singletons_with_partner_all_proteome'] for r in evidence['rows'])
    (out/'evidence-results.tex').write_text(
        f'The reference-evidence annotations contain {reference_singletons} organism--complex singleton pairs across the 50 source proteomes. '
        f'Adding all GO evidence supplies another annotated protein in {same} of these cases on the same reference-annotated proteins, '
        f'and in {whole} cases when all source-proteome proteins are eligible (Supplementary Table~S7). '
        'These are annotation-coverage changes, not validated interaction gains.\n')
    mrows=[]
    for (m,mode),d in memberships.items():
        p,dem=d['promotions'],d['demotions']
        changed_taxa = sum(c['promotions'] + c['demotions'] > 0 for (mm, _, md), c in checks.items() if mm == m and md == mode)
        mrows.append([SHORT[m],mode,changed_taxa,f"{p['total']}/{p['with_curated_reference']}/{p['true_subunit_promotions']}",f"{dem['total']}/{dem['with_curated_reference']}/{dem['curated_subunit_demotions']}",d['unmapped_promotions_with_reference'],d['unmapped_demotions_with_reference']])
    render_table(out,'supp-membership','Curated membership validation. All 50 designated organisms are checked in each row; changed taxa counts those with at least one active-state edit. Triples are total edits / edits with an assessable reference / edits involving a listed subunit. Unmapped counts have an organism--complex reference but unresolved protein accession. Unlisted does not mean biologically incorrect. '+LEGEND,
        ['Model','Mode','Changed taxa','Promotions','Demotions','Unmapped P','Unmapped D'],mrows)
    solver_rows=[]
    for m in MODELS[:3]:
        for mode in ['stage2_standard','stage2_demotion']+[f'soft_stage2_{p}' for p in [.2,.5,1.0]]:
            r=get(m,'cc',mode); counts=r['solver_status_counts']
            if sum(counts.values()) != 50:
                raise ValueError('Missing solver status')
            solver_rows.append([SHORT[m],mode_name(mode),counts.get('OPTIMAL',0),counts.get('FEASIBLE',0),fmt(r['solver_objective_sum']),fmt(r['solver_bound_sum'])])
    render_table(out,'supp-solvers','Bounded Stage 2 optimization. Counts sum to 50 organisms in each row. Objectives and lower bounds are sums across organisms for the indicated objective; differing soft penalties are not the same optimization problem. An optimal restricted model need not be optimal without the candidate restriction. '+LEGEND,
        ['Model','Mode','Optimal','Feasible','Objective sum','Bound sum'],solver_rows)
    budgets = {m: {} for m in MODELS[:3]}
    for (m,a,t), source_record in sources.items():
        if m in budgets and a == 'cc':
            budget = source_record['configuration']['solver_seconds']
            budgets[m][budget] = budgets[m].get(budget,0)+1
    with (out/'supp-solvers.tex').open('a') as fh:
        fh.write('\\par Accepted per-condition search budgets by organism task: '+
            '; '.join(SHORT[m]+': '+', '.join(str(n)+' at '+str(seconds)+' s' for seconds,n in sorted(budgets[m].items())) for m in budgets)+'.\n')
    native=[get(m,a,'stage1') for m in MODELS for a in ASPECTS]
    rates=[r['demotion_active_pct'] for r in native if r['demotion_active_pct'] is not None]
    delta=lambda key:max(abs(get(m,a,'stage1')['metrics'][key]-get(m,a,'before')['metrics'][key]) for m in MODELS for a in ASPECTS)
    burden = ('At the native thresholds, known-taxon Stage~1 demotes '+fmt(min(rates),3)+r'\% to '+fmt(max(rates),3)
        +r'\% of active input pairs across the assessable predictor--aspect conditions (Table~\ref{tab:main-flips}). '
        if rates else r'All active-input denominators are empty at native thresholds; demotion rates are undefined (Table~\ref{tab:main-flips}). ')
    (out/'main-results.tex').write_text(
        burden
        +'The largest absolute change across these conditions is '+fmt(delta('fmax'))
        +r' in $F_{\max}$ and '+fmt(delta('aupr'))+r' in AUPR (Table~\ref{tab:main-metrics}). '
        +'Full threshold-specific counts and metrics are reported in Supplementary Tables~S3 and~S4.\n')
    joint=[r for r in lineage['summary'] if r['aspect']=='joint']
    (out/'lineage-results.tex').write_text('The jointly evaluated GO-supported lineages are exact ancestral prefixes for '+', '.join(SHORT[r['predictor']]+' '+str(r['supported']['exact_prefix_count'])+'/'+str(r['n']) for r in joint)+'. These counts exclude fully unresolved outputs. Supplementary Table~S5 reports all seven ranks, while the accompanying per-organism records distinguish raw assignments, required ranks and cross-aspect compatibility. Published filtering and soft-penalty comparisons are reported in Supplementary Table~S6.\n')
    hard=[get(m,'cc',mode) for m in MODELS[:3] for mode in ['stage2_standard','stage2_demotion']]
    if any(r['known_taxon_consistency'] != 1 or r['complex_coherence'] != 1 for r in hard):
        raise ValueError('Hard Stage 2 structural contract failed')
    nopt=sum(r['solver_status_counts'].get('OPTIMAL',0) for r in hard)
    nfeas=sum(r['solver_status_counts'].get('FEASIBLE',0) for r in hard)
    (out/'complex-results.tex').write_text('Both hard Stage~2 modes satisfy the known-taxon and complex-cardinality checks on all 50 organisms for each of the three entry-date predictors. Of the 300 organism--mode solves, '+str(nopt)+' are proven optimal within their encoded domains and '+str(nfeas)+' are feasible without a proof of optimality at the search limit. Supplementary Tables~S8 and~S9 report curated membership coverage and optimization bounds.\n')
    # The diagram is executable vector source, not a copied rendered figure.
    shutil.copyfile(Path(__file__).parent/'figures/repair-schematic.tex', out/'repair-schematic.tex')
    manifest={'inputs':{str(p):digest(p) for p in required},'generator_sha256':digest(__file__),
              'outputs':{p.name:digest(p) for p in sorted(out.iterdir()) if p.suffix in {'.tex','.tsv'}}}
    (out/'generation-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    destination.mkdir(parents=True, exist_ok=True)
    for path in sorted(out.iterdir(), key=lambda p: p.name == 'generation-manifest.json'):
        path.replace(destination/path.name)
    staging.cleanup()
    print(f'Generated {len(manifest["outputs"])} source artifacts from {len(required)} validated inputs in {destination}')


if __name__ == '__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',required=True,help='Complete result bundle with summaries/, lineage, evidence and membership JSON files')
    ap.add_argument('--output',required=True)
    generate(ap.parse_args())
