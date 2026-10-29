"""Independently validate every serialized hard Stage 2 repair, row by row."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import sys

ROOT = Path(os.environ.get('AUDIT_CODE_DIR', Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT))
from gaef_patches.ontology import Ontology
from taxon_consistency import adjust_ortools as tc
from complex_coherence import adjust_ortools as cc


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def parse(row):
    if not row:
        raise ValueError('Empty prediction row')
    scores = {}
    for token in row[1:]:
        term, value = token.split('|')
        if term in scores:
            raise ValueError(f'Duplicate term {term}')
        score = float(value)
        if not 0 <= score <= 1:
            raise ValueError(f'Invalid score {value}')
        scores[term] = score
    return row[0], scores


def check_pair(before_path, after_path, tau, hierarchy, complexes, forbidden, demotion_only):
    counts, promotions, demotions, proteins = Counter(), 0, 0, 0
    with open(before_path) as f1, open(after_path) as f2:
        for a, b in zip_longest(csv.reader(f1, delimiter='\t'), csv.reader(f2, delimiter='\t')):
            if a is None or b is None:
                raise ValueError('Protein row counts differ')
            protein, before = parse(a)
            protein2, after = parse(b)
            if protein != protein2:
                raise ValueError(f'Protein order or identity changed: {protein}, {protein2}')
            proteins += 1
            active = {g for g, value in after.items() if value > tau}
            if active & forbidden:
                raise ValueError(f'Known-taxon violation: {protein}: {active & forbidden}')
            for child in active:
                missing = hierarchy.get(child, set()) - active
                if missing:
                    raise ValueError(f'Hierarchy violation: {protein}, {child}, {missing}')
            counts.update(active & complexes)
            for term in set(before) | set(after):
                old, new = before.get(term, 0.0), after.get(term, 0.0)
                was, now = old > tau, new > tau
                if was == now:
                    expected = old
                elif now:
                    if demotion_only:
                        raise ValueError(f'Promotion in demotion-only mode: {protein}, {term}')
                    promotions += 1
                    expected = max(old, min(1.0, tau+.001))
                else:
                    demotions += 1
                    expected = min(old, tau)
                if new != expected:
                    raise ValueError(f'Score export mismatch: {protein}, {term}, {old}, {new}, expected {expected}')
    if any(n == 1 for n in counts.values()):
        raise ValueError('Singleton complex remains')
    return {'proteins': proteins, 'promotions': promotions, 'demotions': demotions, 'present_complexes': len(counts)}


def run(args):
    resources, results = Path(args.resources), Path(args.results)
    go = Ontology(resources/'go-basic.obo')
    hierarchy = {g: set(obj['is_a']) for g, obj in go.ont.items()}
    isa_go = Ontology(resources/'go-basic.obo', with_rels=False)
    complexes = cc.get_heteromeric_complexes(go.normalize_terms(cc.load_homodimer_terms(ROOT/'complex_coherence/protein_complexes.tsv')), set(),
        {g: set(obj['is_a']) for g, obj in isa_go.ont.items()})
    only, never = tc.load_constraints(resources/'go_taxon_constraints_extracted_obo.tsv')
    only, never = go.normalize_constraints(only), go.normalize_constraints(never)
    isa, _, unions = tc.load_taxon_hierarchy(resources/'taxon_hierarchy.tsv', resources/'ncbitaxon_hierarchy.tsv')
    taxa = set(Path(args.taxa).read_text().split())
    records=[]
    for model in ['deepgo-se','seq_sim_timeset','sprof-go']:
        tau=.3 if model == 'deepgo-se' else .001
        for taxon in sorted(taxa):
            lineage=tc.get_all_ancestors('NCBITaxon_'+taxon, isa)
            for union, members in unions.items():
                if members & lineage:
                    lineage.add(union)
            forbidden={g for g in set(only)|set(never) if any(t not in lineage for t in only.get(g, ())) or any(t in lineage for t in never.get(g, ()))}
            with gzip.open(results/model/'cc'/f'{taxon}.json.gz','rt') as fh:
                metadata=json.load(fh)
            inputs=list((results/model/'cc'/'stage1').glob(f'optimized_fold_*_taxon_{taxon}.tsv'))
            if len(inputs)!=1:
                raise ValueError(f'Expected one Stage 1 file: {inputs}')
            for mode in ['standard','demotion']:
                after=results/model/'cc'/f'stage2_{mode}'/inputs[0].name
                checked=check_pair(inputs[0],after,tau,hierarchy,complexes,forbidden,mode=='demotion')
                record=next(r for r in metadata['runs'] if r['mode']==f'stage2_{mode}')
                if any(checked[k]!=record[k] for k in ['promotions','demotions']):
                    raise ValueError(f'Exported edits disagree with metadata: {after}')
                if checked['present_complexes'] != record['n_present_complexes']:
                    raise ValueError('Present complex count disagrees with metadata')
                if checked['proteins'] != record['n_proteins']:
                    raise ValueError('Protein count disagrees with metadata')
                solver=record['solver']
                if solver['status'] not in {'OPTIMAL','FEASIBLE'} or solver['bound']>solver['objective']+1e-6:
                    raise ValueError('Invalid optimization status/bounds')
                if solver['optimal'] and abs(solver['bound']-solver['objective'])>1e-6:
                    raise ValueError('Optimal status with nonzero objective gap')
                records.append({'predictor':model,'taxon':taxon,'mode':mode,
                    'stage1_sha256':digest(inputs[0]), 'stage2_sha256':digest(after),
                    'metadata_sha256':digest(results/model/'cc'/f'{taxon}.json.gz'), **checked})
            print(model,taxon,flush=True)
    if len(records)!=300:
        raise ValueError('Expected exactly 300 hard Stage 2 validations')
    Path(args.output).write_text(json.dumps({'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'validated':len(records),'records':records},indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['resources','results','taxa','output']:
        ap.add_argument('--'+name,required=True)
    run(ap.parse_args())
