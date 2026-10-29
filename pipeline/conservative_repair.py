"""Apply evidence-restricted complex repair to a Stage 1 score file.

Evidence JSON files contain objects with protein, term and source fields.
Protection uses the same format and applies to existing active annotations.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from gaef_patches.ontology import Ontology
from complex_coherence import adjust_ortools as cc
from complex_coherence.conservative import repair


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):
            h.update(block)
    return h.hexdigest()


def evidence(path,go):
    result={}
    if path is None:
        return result
    for row in json.loads(path.read_text()):
        term=go.canonical_id(row['term'])
        if term is None:
            raise ValueError(f'Unknown or obsolete evidence term: {row["term"]}')
        if go.ont[term]['namespace'] != 'cellular_component':
            raise ValueError('Complex repair evidence must be in the cellular component aspect')
        pair=term,row['protein']
        if pair in result:
            raise ValueError(f'Duplicate evidence pair: {pair}')
        if not isinstance(row['source'],str) or not row['source'].strip():
            raise ValueError('Evidence requires a nonempty source description')
        result[pair]=row['source']
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['predictions','go','output','audit']:
        ap.add_argument('--'+name,type=Path,required=True)
    for name in ['promotion-evidence','demotion-evidence','protected-annotations']:
        ap.add_argument('--'+name,type=Path)
    ap.add_argument('--complex-exclusions',type=Path,default=ROOT/'complex_coherence/protein_complexes.tsv')
    ap.add_argument('--threshold',type=float,required=True)
    ap.add_argument('--singleton-penalty',type=float,default=1.0)
    ap.add_argument('--time-limit',type=float,default=120)
    args=ap.parse_args()
    input_paths=[p for p in [args.predictions,args.go,args.complex_exclusions,
                 args.promotion_evidence,args.demotion_evidence,args.protected_annotations] if p]
    if args.output.resolve() in {p.resolve() for p in input_paths} or args.audit.resolve() in {p.resolve() for p in input_paths}:
        raise ValueError('Output must not overwrite an input')
    if args.output.resolve()==args.audit.resolve():
        raise ValueError('Prediction output and audit must be separate files')
    go=Ontology(args.go)
    hierarchy={g:set(o['is_a']) for g,o in go.ont.items()}
    isa=Ontology(args.go,with_rels=False)
    exclusions=go.normalize_terms(cc.load_homodimer_terms(args.complex_exclusions))
    complexes=cc.get_heteromeric_complexes(exclusions,set(),{g:set(o['is_a']) for g,o in isa.ont.items()})
    normalization={}
    scores={}
    with args.predictions.open() as f:
        for line in f:
            fields=line.rstrip('\n').split('\t')
            protein=fields[0]
            if not protein or protein in scores:
                raise ValueError('Empty or duplicate protein identifier')
            values={}
            for token in fields[1:]:
                term,value=token.split('|')
                if term in values:
                    raise ValueError('Duplicate protein-term score')
                values[term]=float(value)
            scores[protein]=values
    scores=go.normalize_scores(scores,'cellular_component',normalization)
    scores=go.propagate_scores(scores,normalization)
    final,audit=repair(scores,hierarchy,complexes,evidence(args.promotion_evidence,go),
                       evidence(args.demotion_evidence,go),set(evidence(args.protected_annotations,go)),
                       args.threshold,args.singleton_penalty,args.time_limit)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.audit.parent.mkdir(parents=True,exist_ok=True)
    cc.save_predictions(final,args.output)
    audit.update(normalization=normalization,inputs={str(p):sha(p) for p in input_paths},
                 export_sha256=sha(args.output),source_sha256={str(p.relative_to(ROOT)):sha(p) for p in
                     [Path(__file__).resolve(),ROOT/'complex_coherence/conservative.py',ROOT/'complex_coherence/adjust_ortools.py']})
    args.audit.write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    main()
