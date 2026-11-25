"""Deepest rank that GO taxon constraints can ever force for each organism.

A ranked taxon T can be *required* only if some retained GO term has an
only_in_taxon constraint on T or a descendant of T (or a union all of whose
members are at or below T); never_in constraints only exclude. So the deepest
rank in the true lineage that is an ancestor-or-self of some only_in taxon is an
upper bound on the deepest correct rank, whatever the predictor.
"""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--resources', required=True)
    ap.add_argument('--taxa-file', required=True, help='lines: split taxon')
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology
    from taxon_consistency import adjust_ortools as tc
    from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage, RANKS_7
    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    only, never = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    only = go.normalize_constraints(only)
    isa, dis, uni = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    parents, ranks, names = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
    only_taxa = Counter(t for ts in only.values() for t in ts)

    def anc(t):
        out, cur = set(), t
        while cur and cur not in out:
            out.add(cur)
            cur = parents.get(cur)
        return out
    rows = []
    for line in open(a.taxa_file):
        split, taxon = line.split()
        lineage = get_true_lineage(taxon, parents)
        lset = set(lineage)
        forcing = set()  # lineage nodes that some only_in constraint (incl. unions) can force
        for t, n in only_taxa.items():
            if t.startswith('NCBITaxon_Union_'):
                members = uni.get(t, set())
                if not members:
                    continue
                common = set.intersection(*(anc(m) for m in members))
                forcing |= common & lset
            else:
                forcing |= anc(t) & lset
        ceiling = 'none'
        for r in RANKS_7:
            true_r = next((x for x in lineage if ranks.get(x) == r), None)
            if true_r and true_r in forcing:
                ceiling = r
        dom = next((names.get(x, x) for x in lineage if ranks.get(x) == 'domain'), None)
        rows.append({'split': split, 'taxon': taxon, 'domain': dom, 'ceiling': ceiling,
                     'forcing_ranked': sorted((ranks.get(x), names.get(x, x)) for x in forcing if ranks.get(x) in RANKS_7)})
    summary = {}
    for split in sorted({r['split'] for r in rows}):
        summary[split] = {'ceiling': dict(Counter(r['ceiling'] for r in rows if r['split'] == split)),
                          'domain': dict(Counter(r['domain'] for r in rows if r['split'] == split))}
    Path(a.output).write_text(json.dumps({'summary': summary, 'rows': rows}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
