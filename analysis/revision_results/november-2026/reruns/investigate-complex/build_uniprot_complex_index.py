"""Per-organism UniProt GO complex-membership index (silver reference).

For each organism DAT (paths taken from evidence-sensitivity.json), record for
each entry name the operational heteromeric complex terms it is annotated to
(propagated over is_a+part_of), with the strongest evidence class:
  curated (EXP/IDA/IPI/IMP/IGI/IEP/HTP/HDA/HMP/HGI/HEP/TAS/IC) > noniea (other
  non-IEA, e.g. ISS/ISO/IBA) > iea.
This is NOT independent truth (IEA/ISS may share sources with predictors); it
is reported separately from Complex Portal, by evidence class.
"""
import argparse, gzip, json, sys
from pathlib import Path

DIRECT = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
CURATED = DIRECT | {'TAS', 'IC'}
RANK = {'curated': 3, 'noniea': 2, 'iea': 1}

def entries(path):
    entry, reviewed, annotations = None, None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('ID   '):
                f = line.split(); entry = f[1]; reviewed = f[2].startswith('Reviewed')
            elif line.startswith('DR   GO;'):
                fields = [x.strip() for x in line[5:].split(';')]
                if len(fields) >= 4:
                    annotations.append((fields[1], fields[3].split(':')[0]))
            elif line.startswith('//'):
                yield entry, reviewed, annotations
                entry, reviewed, annotations = None, None, []

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True); ap.add_argument('--go', required=True)
    ap.add_argument('--evidence', required=True); ap.add_argument('--out', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology
    from complex_coherence.adjust_ortools import load_homodimer_terms, get_heteromeric_complexes
    go = Ontology(a.go)
    isa = Ontology(a.go, with_rels=False)
    complexes = get_heteromeric_complexes(go.normalize_terms(load_homodimer_terms(Path(a.code)/'complex_coherence/protein_complexes.tsv')), set(),
                                          {g: set(o['is_a']) for g, o in isa.ont.items()})
    rows = json.load(open(a.evidence))['rows']
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    cache = {}
    for r in rows:
        idx, rev = {}, {}
        for entry, reviewed, ann in entries(r['file']):
            rev[entry] = reviewed
            d = {}
            for term, ev in ann:
                t = go.canonical_id(term) or term
                if t not in cache:
                    cache[t] = (go.get_ancestors(t) | {t}) & complexes if go.has_term(t) else set()
                cls = 'curated' if ev in CURATED else ('iea' if ev == 'IEA' else 'noniea')
                for g in cache[t]:
                    if RANK[cls] > RANK.get(d.get(g), 0):
                        d[g] = cls
            if d:
                idx[entry] = d
        with gzip.open(out/f"{r['taxon']}.json.gz", 'wt') as fh:
            json.dump({'taxon': r['taxon'], 'file': r['file'], 'n_entries': len(rev), 'reviewed': sorted(e for e, v in rev.items() if v), 'index': idx}, fh)
        print(r['taxon'], len(rev), len(idx), flush=True)

if __name__ == '__main__':
    main()
