"""Closure-based taxon audit of curated reference annotations (corrected pipeline).

Uses the same ontology handling (primary IDs, is_a + part_of closure) and the
same known-lineage rule as analysis/revision_audit_run.py: the lineage contains
the organism, its NCBI ancestors and every taxon union with a member in the
lineage. An organism is inconsistent if the GO closure of its annotations
contains a term whose only_in_taxon set misses the lineage or whose
never_in_taxon set hits it. The bundled GAEF pairwise screen is reported
separately per organism.
"""
import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--resources', required=True)
    ap.add_argument('--annotations', required=True)
    ap.add_argument('--taxa', help='file with one taxon id per line (restrict)')
    ap.add_argument('--output', required=True)
    ap.add_argument('--top', type=int, default=20)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology
    from taxon_consistency import adjust_ortools as tc
    from pipeline.genome_scores_evaluator import load_annotations, check_consistency_fast

    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    in_tax, never_tax = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    in_tax, never_tax = go.normalize_constraints(in_tax), go.normalize_constraints(never_tax)
    tax_isa, _, tax_union = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    keep = set(Path(a.taxa).read_text().split()) if a.taxa else None
    constraints_tsv = Path(a.code) / 'data/constraints/taxon_constraints.tsv'

    anc_cache = {}
    per_org = {}
    pair_orgs, pair_ann = defaultdict(set), defaultdict(int)
    for f in sorted(Path(a.annotations).glob('annots_taxon_*.tsv')):
        org = re.search(r'annots_taxon_(\d+)\.tsv$', f.name).group(1)
        if keep is not None and org not in keep:
            continue
        ann = {p: go.normalize_terms(ts) for p, ts in load_annotations(f).items()}
        ann = {p: ts for p, ts in ann.items() if ts}
        lineage = {tc.normalize_taxon_id(org)}
        lineage.update(tc.get_all_ancestors_of_list(tax_isa, list(lineage)))
        for u, members in tax_union.items():
            if members & lineage:
                lineage.add(u)
        def bad(g):
            only = in_tax.get(g, ())
            if only and not any(t in lineage for t in only):
                return ('only_in_taxon', ','.join(sorted(only)))
            hit = [t for t in never_tax.get(g, ()) if t in lineage]
            return ('never_in_taxon', ','.join(sorted(hit))) if hit else None
        genome = set()
        violated = defaultdict(int)
        for p, ts in ann.items():
            closed = set()
            for t in ts:
                if t not in anc_cache:
                    anc_cache[t] = go.get_ancestors(t)
                closed |= anc_cache[t]
            genome |= closed
            for g in closed:
                b = bad(g)
                if b:
                    violated[(g,) + b] += 1
        for key, n in violated.items():
            pair_orgs[key].add(org)
            pair_ann[key] += n
        per_org[org] = {'n_proteins': len(ann), 'n_terms': len(genome),
                        'closure_violation': bool(violated),
                        'violating_terms': sorted({k[0] for k in violated}),
                        'gaef_pairwise_consistent': bool(check_consistency_fast({'genome': genome}, constraints_tsv, ontology=go)) if genome else None}
    ranked = sorted(pair_orgs, key=lambda k: (len(pair_orgs[k]), pair_ann[k]), reverse=True)
    names = {}
    for g, _, _ in ranked[:a.top]:
        try:
            names[g] = go.get_term(g)['name']
        except Exception:
            names[g] = ''
    out = {'per_organism': per_org,
           'n_organisms': len(per_org),
           'n_closure_violating': sum(v['closure_violation'] for v in per_org.values()),
           'n_gaef_inconsistent': sum(v['gaef_pairwise_consistent'] is False for v in per_org.values()),
           'top_pairs': [{'term': k[0], 'name': names.get(k[0], ''), 'constraint': k[1], 'taxa': k[2],
                          'n_organisms': len(pair_orgs[k]), 'n_annotations': pair_ann[k]} for k in ranked[:a.top]]}
    Path(a.output).write_text(json.dumps(out, indent=1))
    print(json.dumps({k: v for k, v in out.items() if k != 'per_organism'}, indent=1))


if __name__ == '__main__':
    main()
