"""Summarize variant outputs: membership support of edits, coherence, test-protein effects.

Membership references:
  CP   Complex Portal index (as validate_stage2.py: assessable = organism-complex
       key indexed and accession mapped; hit = listed subunit)
  UNI  UniProt GO annotation of the same protein to the same complex term
       (propagated), any evidence / non-IEA / curated. Absence is not a negative,
       so UNI-precision is a lower bound and only comparable across policies.
Edits are complex-term pairs (g in operational complexes) in the Stage1->variant diff.
Test-protein effects: promoted/demoted pairs on temporal test proteins, checked
against the propagated reference annotations stored in the accepted result JSON.
"""
import argparse, gzip, json
from collections import defaultdict, Counter
from pathlib import Path
import common_ic as C


def main():
    ap = argparse.ArgumentParser()
    for k in ['code', 'go', 'accepted', 'cp_index', 'uniprot_dir', 'idmap', 'vdir', 'out']:
        ap.add_argument('--'+k, required=True)
    ap.add_argument('--taxa', required=True)
    a = ap.parse_args()
    go, H, isa_h, complexes, cc = C.resources(a.code, a.go)
    cp = C.load_cp_index(a.cp_index, go, complexes)
    mapping = json.load(open(a.idmap))
    taxa = a.taxa.split(',')
    anc_cache = {}
    def ancs(g):
        if g not in anc_cache:
            anc_cache[g] = go.get_ancestors(g) | {g}
        return anc_cache[g]
    summary = {}
    for pred in C.PREDS:
        agg = defaultdict(Counter)
        per_case = []
        for T in taxa:
            f = Path(a.vdir)/pred/f'{T}.json.gz'
            if not f.exists():
                agg['_missing']['taxa'] += 1
                continue
            d = json.load(gzip.open(f, 'rt'))
            uni = C.load_uniprot_index(Path(a.uniprot_dir)/f'{T}.json.gz')['index']
            with gzip.open(Path(a.accepted)/pred/'cc'/f'{T}.json.gz', 'rt') as fh:
                res = json.load(fh)
            ref = {p: set().union(*(ancs(g) for g in ts)) if ts else set() for p, ts in res['annotations'].items()}
            del res
            for v, r in d['variants'].items():
                A = agg[v]
                A['taxa'] += 1
                A['initial_singletons'] += d['n_singletons']
                A['present_after'] += r['n_present']
                A['singletons_after'] += r['n_singletons_after']
                A.update('outcome_'+o for o in r['outcome'].values())
                A['optimal'] += r['solver'].get('optimal', True) is True
                for kind, pairs in [('prom', r['promotions']), ('dem', r['demotions'])]:
                    A[kind+'_all'] += len(pairs)
                    for g, p in pairs:
                        if p in ref:
                            A[kind+'_test'] += 1
                            A[kind+'_test_inref'] += g in ref[p]
                        if g not in complexes:
                            continue
                        A[kind+'_cx'] += 1
                        acc = C.protein_to_accession(p, mapping)
                        key = f'{T}::{g}'
                        if key in cp:
                            if acc is None:
                                A[kind+'_cp_unmapped'] += 1
                            else:
                                A[kind+'_cp_assess'] += 1
                                A[kind+'_cp_listed'] += acc in cp[key]
                        cls = uni.get(p, {}).get(g)
                        A[kind+'_uni_any'] += cls is not None
                        A[kind+'_uni_noniea'] += cls in ('noniea', 'curated')
                        A[kind+'_uni_curated'] += cls == 'curated'
                # singleton-level: did the repair of X add a UniProt/CP member?
                prom_by_X = defaultdict(set)
                for g, p in r['promotions']:
                    prom_by_X[g].add(p)
                for X, o in r['outcome'].items():
                    if o == 'promoted':
                        A['X_prom'] += 1
                        A['X_prom_uni'] += any(X in uni.get(p, {}) for p in prom_by_X[X])
                        key = f'{T}::{X}'
                        if key in cp:
                            accs = {C.protein_to_accession(p, mapping) for p in prom_by_X[X]} - {None}
                            if accs:
                                A['X_prom_cp_assess'] += 1
                                A['X_prom_cp_listed'] += bool(accs & cp[key])
                    if o == 'demoted':
                        L = d['lone'][X]
                        A['X_dem'] += 1
                        A['X_dem_lone_uni'] += X in uni.get(L, {})
                        key = f'{T}::{X}'
                        acc = C.protein_to_accession(L, mapping)
                        if key in cp and acc:
                            A['X_dem_cp_assess'] += 1
                            A['X_dem_cp_listed'] += acc in cp[key]
            # evidence-level precision (independent of solver)
            E = agg['_evidence']
            for X, evs in d['evidence'].items():
                key = f'{T}::{X}'
                for p, types in evs.items():
                    for t in types + ['any']:
                        E[t+'_pairs'] += 1
                        E[t+'_uni'] += X in uni.get(p, {})
                        acc = C.protein_to_accession(p, mapping)
                        if key in cp and acc:
                            E[t+'_cp_assess'] += 1
                            E[t+'_cp_listed'] += acc in cp[key]
                E['singletons'] += 1
                E['singletons_with_evidence'] += bool(evs)
        summary[pred] = {k: dict(v) for k, v in agg.items()}
    Path(a.out).write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
