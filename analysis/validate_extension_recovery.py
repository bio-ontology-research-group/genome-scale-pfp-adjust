"""Reconstruct masked-recovery eligibility and independently verify outcomes."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from analysis.revision_extensions import resources, load_selected, sha
from analysis.complex_validation.validate_stage2 import protein_to_accession

ACC=re.compile(r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$')
MODES={'standard','demotion_only','conservative','conservative_corroborated','highest_score','highest_score_corroborated'}
PREDICTORS={'deepgo-se','seq_sim_timeset','sprof-go'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['output-root','results','resources','id-map','complexportal','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    go,_,complexes=resources(args.resources)
    mapping=json.loads(args.id_map.read_text())
    expected={(p,t) for p,a,t in (line.split() for line in (ROOT/'reproducibility/revision_tasks.txt').read_text().splitlines()) if p in PREDICTORS}
    paths={(p.parent.name,p.stem):p for p in args.output_root.glob('recovery/*/*.json')}
    assert set(paths)==expected
    records=[];inputs={};case_count=0
    for taxon in sorted({t for p,t in expected}):
        cp=args.complexportal/f'{taxon}.tsv'
        results={p:json.loads(paths[p,taxon].read_text()) for p in PREDICTORS}
        if not cp.exists():
            assert all(not d['records'] and d['exclusion']=='no_complexportal_taxon' for d in results.values())
            records.extend(dict(predictor=p,taxon=taxon,cases=0,result_sha256=sha(paths[p,taxon]),status='PASS') for p in sorted(PREDICTORS))
            continue
        candidates=set()
        with cp.open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                if row['Taxonomy identifier'].split('(')[0].strip()!=taxon:
                    continue
                tokens=[t.split('(')[0].strip() for t in row['Identifiers (and stoichiometry) of molecules in complex'].split('|')]
                accessions=tuple(sorted({t for t in tokens if ACC.fullmatch(t)}))
                if len(accessions)==2:
                    for term in set(re.findall(r'GO:\d{7}',row['Go Annotations'])):
                        canonical=go.canonical_id(term)
                        if canonical in complexes:
                            candidates.add((canonical,accessions))
        scores={};source_paths={}
        for p in PREDICTORS:
            files=list((args.results/p/'cc'/'stage1').glob(f'*_taxon_{taxon}.tsv'))
            assert len(files)==1
            source_paths[p]=files[0]
            scores[p]=load_selected(files[0],{g for g,a in candidates})
            inputs[str(files[0])]=sha(files[0])
        inputs[str(cp)]=sha(cp)
        for p,data in sorted(results.items()):
            original=scores[p];tau=.3 if p=='deepgo-se' else .001
            assert data['threshold']==tau and data['n_proteins']==len(original)
            assert data['inputs'][str(source_paths[p])]==inputs[str(source_paths[p])]
            assert data['inputs'][str(cp)]==inputs[str(cp)]
            assert data['inputs'][str(args.id_map)]==sha(args.id_map)
            for other in PREDICTORS-{p}:
                assert data['supporting_inputs'][str(source_paths[other])]==inputs[str(source_paths[other])]
            inverse=defaultdict(list)
            for protein in original:
                acc=protein_to_accession(protein,mapping)
                if acc: inverse[acc].append(protein)
            cases=set()
            for term,accessions in candidates:
                if all(len(inverse[a])==1 for a in accessions):
                    a,b=(inverse[x][0] for x in accessions)
                    cases.update([(term,a,b),(term,b,a)])
            actual={(r['term'],r['seed'],r['hidden']) for r in data['records']}
            assert cases==actual and len(cases)==len(data['records'])
            for r in data['records']:
                term,seed,hidden=r['term'],r['seed'],r['hidden']
                primary={protein for protein,s in original.items() if protein!=seed and s.get(term,0)>tau}
                corroborated={protein for protein in original if protein!=seed and all(
                    scores[other].get(protein,{}).get(term,0)>(.3 if other=='deepgo-se' else .001)
                    for other in PREDICTORS-{p})}
                assert r['original_hidden_score']==original[hidden].get(term,0)
                assert r['hidden_has_evidence']==(hidden in primary)
                assert r['hidden_has_corroboration']==(hidden in corroborated)
                assert r['eligible_candidates']==len(primary)
                masked={protein:max(.9,tau+.001) if protein==seed else tau*s.get(term,0) for protein,s in original.items()}
                assert {protein for protein,s in masked.items() if s>tau}=={seed}
                top50=set(sorted(original,key=lambda protein:(-masked[protein],protein))[:50])
                assert set(r['outputs'])==MODES
                for mode,out in r['outputs'].items():
                    members=set(out['members'])
                    assert len(members)==len(out['members']) and members<=set(original)
                    assert out['hidden_recovered']==(hidden in members)
                    assert out['seed_removed']==(seed not in members)
                    assert out['unresolved']==(len(members)==1)
                    assert set(out['other_proposals'])==members-{seed,hidden}
                    if mode=='standard':
                        assert len(members)!=1 and members<=top50
                    elif mode=='demotion_only':
                        assert not members
                    else:
                        allowed=corroborated if mode.endswith('corroborated') else primary
                        assert seed in members and members-{seed}<=allowed
                        if mode.startswith('highest_score'):
                            ranked=sorted(allowed,key=lambda protein:(-original[protein].get(term,0),protein))
                            assert members=={seed}|set(ranked[:1])
                    if not mode.startswith('highest_score'):
                        assert out['solver']['status'] in ['OPTIMAL','FEASIBLE']
            case_count+=len(cases)
            records.append(dict(predictor=p,taxon=taxon,cases=len(cases),result_sha256=sha(paths[p,taxon]),status='PASS'))
    inputs[str(args.id_map)]=sha(args.id_map)
    args.output.write_text(json.dumps(dict(status='PASS',n_tasks=len(records),n_cases=case_count,records=records,
                                          inputs=inputs,validator_sha256=sha(Path(__file__))),indent=2)+'\n')
    print(f'PASS: {len(records)} tasks, {case_count} recovery cases')

if __name__=='__main__':
    main()
