"""Independent streaming checks of conservative exports and curated coverage."""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.revision_extensions import resources, sha, load_selected, corroborating_evidence
from analysis.complex_validation.validate_stage2 import protein_to_accession
from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage
from taxon_consistency import adjust_ortools as tc


def row(line):
    cols = line.rstrip('\n').split('\t')
    values = {}
    for token in cols[1:]:
        if '|' in token:
            g,v = token.split('|')
            if g in values:
                raise ValueError('Duplicate score pair')
            values[g] = float(v)
    return cols[0], values


def validate(args, result, go, hierarchy, complexes, only, never, isa, unions, mapping, reference):
    audit = json.loads(result.read_text())
    predictor, taxon = audit['predictor'], audit['taxon']
    tau = audit['threshold']
    inputs = sorted((args.results/predictor/'cc'/'stage1').glob(f'*_taxon_{taxon}.tsv'))
    if len(inputs)!=1:
        raise ValueError('Missing input')
    original, exported = inputs[0], result.with_suffix('.tsv')
    assert sha(original)==audit['input_sha256'] and sha(exported)==audit['export_sha256']
    lineage = {tc.normalize_taxon_id(taxon)} | tc.get_all_ancestors_of_list(isa, [tc.normalize_taxon_id(taxon)])
    lineage |= {u for u, members in unions.items() if members & lineage}
    before_counts, after_counts = Counter(), Counter()
    changed, proteins = set(), set()
    stats = Counter()
    with open(original) as left, open(exported) as right:
        for a,b in itertools.zip_longest(left,right):
            assert a is not None and b is not None, 'Protein count differs'
            p, x = row(a); q,y = row(b)
            assert p==q and p not in proteins, 'Protein identity/order differs'
            proteins.add(p)
            before = {g for g,v in x.items() if v>tau}
            after = {g for g,v in y.items() if v>tau}
            assert before <= after, 'Conservative mode deleted an active annotation'
            assert all(parent in after for g in after for parent in hierarchy.get(g,())), 'Hierarchy violation'
            assert all(set(only.get(g,())) <= lineage and not set(never.get(g,())) & lineage for g in after), 'Taxon violation'
            before_counts.update(before & complexes); after_counts.update(after & complexes)
            for g in x.keys() | y.keys():
                if x.get(g,0)==y.get(g,0):
                    continue
                assert g in after-before, 'Unrelated score changed'
                assert y[g]==max(x.get(g,0), min(1,tau+.001)), 'Incorrect score reconstruction'
                changed.add((g,p))
                if g in complexes:
                    stats['complex_promotions'] += 1
                    accession = protein_to_accession(p, mapping)
                    listed = reference.get(f'{taxon}::{g}')
                    if accession and listed is not None:
                        stats['assessable_promotions'] += 1
                        stats['listed_promotions'] += int(accession in listed)
                    elif accession is None:
                        stats['unresolved_accessions'] += 1
    assert changed=={tuple(pair) for pair in audit['promotions']}
    assert not audit['demotions']
    singletons = {g for g,n in before_counts.items() if n==1}
    if audit['evidence_policy'].startswith('agreement'):
        evidence_args = argparse.Namespace(results=args.results, predictor=predictor, taxon=taxon)
        evidence, sources = corroborating_evidence(evidence_args, singletons)
        assert sources == audit['supporting_inputs']
        allowed = {(a,p) for g,p in evidence for a in go.get_ancestors(g)}
    else:
        selected = load_selected(original, singletons)
        allowed = {(a,p) for p,terms in selected.items() for g,v in terms.items()
                   if 0<v<=tau for a in go.get_ancestors(g)}
    assert changed <= allowed, 'Promotion lacks authorized evidence or its ancestor implication'
    unresolved = {g for g,n in after_counts.items() if n==1}
    assert unresolved=={r['term'] for r in audit['unresolved']}
    assert len(singletons)==audit['initial_singletons']
    assert unresolved <= singletons, 'Promotion-only mode created a singleton'
    return dict(predictor=predictor,taxon=taxon,mode=result.parents[1].name,
                proteins=len(proteins), promotions=len(changed), initial_singletons=len(singletons),
                unresolved=len(unresolved), membership=dict(stats), status='PASS',
                input_sha256=sha(original), export_sha256=sha(exported), result_sha256=sha(result))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ['output-root','results','resources','id-map','subunit-index','output']:
        ap.add_argument('--'+key,type=Path,required=True)
    ap.add_argument('--worker',type=int,default=0)
    ap.add_argument('--workers',type=int,default=1)
    args=ap.parse_args()
    go,hierarchy,complexes=resources(args.resources)
    only,never=tc.load_constraints(args.resources/'go_taxon_constraints_extracted_obo.tsv')
    only,never=go.normalize_constraints(only),go.normalize_constraints(never)
    isa,_,unions=tc.load_taxon_hierarchy(args.resources/'taxon_hierarchy.tsv',args.resources/'ncbitaxon_hierarchy.tsv')
    mapping=json.loads(args.id_map.read_text())
    raw=json.loads(args.subunit_index.read_text())
    reference={}
    for key,values in raw.items():
        taxon,term=key.split('::'); term=go.canonical_id(term)
        if term in complexes:
            reference.setdefault(f'{taxon}::{term}',set()).update(values)
    paths=sorted(list(args.output_root.glob('conservative/*/*.json'))+list(args.output_root.glob('corroborated/*/*.json')))
    if len(paths)!=300:
        raise ValueError(f'Expected 300 exports, found {len(paths)}')
    records=[validate(args,p,go,hierarchy,complexes,only,never,isa,unions,mapping,reference)
             for i,p in enumerate(paths) if i%args.workers==args.worker]
    provenance = dict(index_sha256=sha(args.subunit_index), mapping_sha256=sha(args.id_map),
                      go_sha256=sha(args.resources/'go-basic.obo'),
                      complex_exclusions_sha256=sha(ROOT/'complex_coherence/protein_complexes.tsv'))
    args.output.write_text(json.dumps(dict(records=records, provenance=provenance,
                                        validator_sha256=sha(Path(__file__))),indent=2)+'\n')

if __name__=='__main__':
    main()
