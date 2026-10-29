"""Aggregate only a complete, independently checked extension experiment."""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.revision_extensions import sha
from analysis.evaluate_lineage_inference import RANKS_7


def write_tsv(path, rows):
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['lineage-root','complex-root','validation-root','recovery-validation','baseline-root','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    expected={(p,t) for p,a,t in (line.split() for line in (ROOT/'reproducibility/revision_tasks.txt').read_text().splitlines())}
    entry={(p,t) for p,t in expected if p in ['deepgo-se','seq_sim_timeset','sprof-go']}
    assert len(expected)==588 and len(entry)==150
    inputs={}
    def load(path):
        inputs[str(path)]=sha(path)
        return json.loads(path.read_text())
    qc={}
    qc_provenance=None
    for path in sorted(args.validation_root.glob('exports-*.json')):
        validation=load(path)
        if qc_provenance is None:
            qc_provenance=validation['provenance']
        assert validation['provenance']==qc_provenance
        for r in validation['records']:
            key=(r['mode'],r['predictor'],r['taxon'])
            assert key not in qc and r['status']=='PASS'
            qc[key]=r
    assert set(qc)=={(m,p,t) for m in ['conservative','corroborated'] for p,t in entry}, 'Incomplete export validation'
    recovery_validation=load(args.recovery_validation)
    assert recovery_validation['status']=='PASS' and recovery_validation['n_tasks']==150
    recovery_qc={(r['predictor'],r['taxon']):r for r in recovery_validation['records']}
    assert set(recovery_qc)==entry
    lineage_groups=defaultdict(list); rank_rows=[]
    for p,t in sorted(expected):
        data=load(args.lineage_root/'lineage'/p/f'{t}.json')
        assert (data['predictor'],data['taxon'])==(p,t)
        assert {r['aspect'] for r in data['records']}=={'cc','mf','bp','joint'}
        for r in data['records']:
            audit=r['uncertainty']
            assert audit['near_optimal_units']==0
            assert set(audit['required']) <= set(r['conditional_required'])
            assert all(audit['checks'][taxon]=='INFEASIBLE' for taxon in audit['required'])
            lineage_groups[p,r['aspect']].append(r)
            conditional,robust=r['conditional_evaluation'],r['robust_evaluation']
            assert conditional['true_ranks']==robust['true_ranks']
            for rank in RANKS_7:
                truth=robust['true_ranks'].get(rank,'')
                c=conditional['inferred_ranks'].get(rank,[]);v=robust['inferred_ranks'].get(rank,[])
                rank_rows.append(dict(predictor=p,taxon=t,aspect=r['aspect'],rank=rank,true_taxon=truth,
                                      conditional=';'.join(c),robust=';'.join(v),
                                      conditional_resolved=bool(c),robust_resolved=bool(v),
                                      conditional_correct=bool(truth) and c==[truth],
                                      robust_correct=bool(truth) and v==[truth]))
    lineage=[]
    for (p,aspect),records in sorted(lineage_groups.items()):
        lineage.append(dict(predictor=p,aspect=aspect,n=len(records),
            optimal=sum(r['uncertainty']['status']=='OPTIMAL' for r in records),
            queries=sum(len(r['uncertainty'].get('checks',{})) for r in records),
            undetermined=sum(len(r['uncertainty']['undetermined']) for r in records),
            conditional_claims=sum(len(r['conditional_required']) for r in records),
            robust_claims=sum(len(r['uncertainty']['required']) for r in records),
            cases_losing_claims=sum(bool(r['removed_conditional']) for r in records),
            conditional_correct_prefix=sum(r['conditional_evaluation']['is_exact_prefix'] for r in records),
            robust_correct_prefix=sum(r['robust_evaluation']['is_exact_prefix'] for r in records),
            robust_unresolved=sum(not r['uncertainty']['required'] for r in records),
            robust_wrong=sum(r['robust_evaluation']['first_divergence_rank']!='none' for r in records)))
    complex_groups=defaultdict(list); recovery_groups=defaultdict(list); case_rows=[]; complex_rows=[]
    recovery_exclusions=Counter(); recovery_coverage=[]
    for p,t in sorted(entry):
        for mode in ['conservative','corroborated']:
            path=args.complex_root/mode/p/f'{t}.json'
            data=load(path); q=qc[mode,p,t]
            assert inputs[str(path)]==q['result_sha256']
            assert data['export_sha256']==q['export_sha256'] and data['input_sha256']==q['input_sha256']
            assert len(data['promotions'])==q['promotions'] and not data['demotions']
            complex_groups[p,mode].append((data,q))
            complex_rows.append(dict(predictor=p,taxon=t,mode=mode,initial_singletons=q['initial_singletons'],
                unresolved=q['unresolved'],promotions=q['promotions'],
                assessable=q['membership'].get('assessable_promotions',0),listed=q['membership'].get('listed_promotions',0),
                solver_status=data['solver']['status'],objective=data['solver']['objective'],bound=data['solver']['bound']))
        recovery_path=args.complex_root/'recovery'/p/f'{t}.json'
        data=load(recovery_path)
        assert inputs[str(recovery_path)]==recovery_qc[p,t]['result_sha256']
        assert len(data['records'])==recovery_qc[p,t]['cases']
        recovery_coverage.append(dict(predictor=p,taxon=t,cases=len(data['records']),
                                      exclusion=data.get('exclusion',''),details=data.get('excluded',{})))
        recovery_exclusions.update(data.get('excluded',{}))
        for r in data['records']:
            assert set(r['outputs'])=={'standard','demotion_only','conservative','conservative_corroborated','highest_score','highest_score_corroborated'}
            assert r['seed']!=r['hidden']
            for mode,output in r['outputs'].items():
                assert output['hidden_recovered']==(r['hidden'] in output['members'])
                assert output['seed_removed']==(r['seed'] not in output['members'])
                assert output['unresolved']==(len(output['members'])==1)
                if mode.startswith('conservative') or mode.startswith('highest_score'):
                    assert not output['seed_removed']
                recovery_groups[p,mode].append((t,r,output))
                case_rows.append(dict(predictor=p,taxon=t,term=r['term'],seed=r['seed'],hidden=r['hidden'],mode=mode,
                    original_hidden_score=r['original_hidden_score'],hidden_has_evidence=r['hidden_has_evidence'],
                    hidden_has_corroboration=r['hidden_has_corroboration'],
                    recovered=output['hidden_recovered'],seed_removed=output['seed_removed'],unresolved=output['unresolved'],
                    other_proposals=';'.join(output['other_proposals']),members=';'.join(output['members']),
                    solver_status=output.get('solver',{}).get('status','not_applicable')))
    complexes=[]
    for (p,mode),records in sorted(complex_groups.items()):
        complexes.append(dict(predictor=p,mode=mode,n=len(records),
            initial_singletons=sum(q['initial_singletons'] for d,q in records),
            unresolved=sum(q['unresolved'] for d,q in records),
            promotions=sum(q['promotions'] for d,q in records),demotions=0,
            assessable=sum(q['membership'].get('assessable_promotions',0) for d,q in records),
            listed=sum(q['membership'].get('listed_promotions',0) for d,q in records),
            listed_demotions=0,optimal=sum(d['solver']['optimal'] for d,q in records),
            objective=sum(d['solver']['objective'] for d,q in records),bound=sum(d['solver']['bound'] for d,q in records)))
    # Existing hard results remain independently validated comparators, not reruns.
    for p in sorted({p for p,t in entry}):
        baseline=load(args.baseline_root/'summaries'/f'{p}_cc.json')
        for mode in ['standard','demotion']:
            label='stage2_standard' if mode=='standard' else 'stage2_demotion'
            rows=[r for r in baseline['rows'] if r['mode']==label]
            assert len(rows)==1, label
            r=rows[0]; m=load(args.baseline_root/f'membership_{p}_{mode}.json')
            assert all(m['provenance'][k]==v for k,v in qc_provenance.items())
            baseline_inputs={pair['taxon']:pair['stage1'] for pair in m['provenance']['pairs']}
            assert all(qc[policy,p,t]['input_sha256']==baseline_inputs[t]
                       for policy in ['conservative','corroborated'] for predictor,t in entry if predictor==p)
            initial=next(x['initial_singletons'] for x in complexes if x['predictor']==p)
            complexes.append(dict(predictor=p,mode='hard_'+mode,n=r['n_organisms'],initial_singletons=initial,unresolved=0,
                promotions=r['promotions'],demotions=r['demotions'],assessable=m['promotions']['with_curated_reference'],
                listed=m['promotions']['true_subunit_promotions'],listed_demotions=m['demotions']['curated_subunit_demotions'],
                optimal=r['solver_status_counts'].get('OPTIMAL',0),objective=r['solver_objective_sum'],bound=r['solver_bound_sum']))
    recovery=[]
    for (p,mode),records in sorted(recovery_groups.items()):
        recovery.append(dict(predictor=p,mode=mode,cases=len(records),taxa=len({t for t,r,o in records}),
            hidden_supported=sum(r['hidden_has_evidence'] for t,r,o in records),
            hidden_corroborated=sum(r['hidden_has_corroboration'] for t,r,o in records),
            recovered=sum(o['hidden_recovered'] for t,r,o in records),
            seed_removed=sum(o['seed_removed'] for t,r,o in records),
            unresolved=sum(o['unresolved'] for t,r,o in records),
            other_proposals=sum(len(o['other_proposals']) for t,r,o in records),
            optimal=sum(o.get('solver',{}).get('status')=='OPTIMAL' for t,r,o in records)))
    args.output.mkdir(parents=True,exist_ok=True)
    for name,rows in [('lineage-summary',lineage),('lineage-ranks',rank_rows),('complex-summary',complexes),
                      ('complex-cases',complex_rows),('recovery-summary',recovery),('recovery-cases',case_rows)]:
        write_tsv(args.output/(name+'.tsv'),rows)
    report=dict(lineage=lineage,complexes=complexes,recovery=recovery,recovery_coverage=recovery_coverage,
                recovery_exclusions=dict(recovery_exclusions),inputs=inputs,
                counts=dict(lineage_tasks=len(expected),lineage_outcomes=sum(x['n'] for x in lineage),
                            rank_rows=len(rank_rows),validated_exports=len(qc),recovery_tasks=len(entry)),
                source_sha256=sha(Path(__file__)))
    (args.output/'summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report['counts']))

if __name__=='__main__':
    main()
