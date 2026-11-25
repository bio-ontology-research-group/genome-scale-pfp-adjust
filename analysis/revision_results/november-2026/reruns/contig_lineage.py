"""Joint taxon recovery from sampled contigs with the corrected Stage 1 solver.

For each contig length, sample contigs at random genome positions (seed 42, as
in the original analysis), take the proteins whose CDS overlaps the contig,
and run Stage 1 without a taxon on their canonical, ancestor-propagated scores.
A ranked taxon (domain..species) counts as inferred only if every optimal GO
repair requires it (feasibility check of its negation at the proven optimum),
as in the revised lineage evaluation. Each contig is classified as
no_annotations, unresolved (no required ranked taxon), wrong (a required taxon
outside the true lineage) or correct to its deepest correct rank.
"""
import argparse
import csv
import json
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--resources', required=True)
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--cds', required=True)
    ap.add_argument('--genome-size', type=int, required=True)
    ap.add_argument('--taxon', required=True)
    ap.add_argument('--aspect', required=True)
    ap.add_argument('--tau', type=float, required=True)
    ap.add_argument('--lengths', default='1000,10000,100000,1000000')
    ap.add_argument('--samples', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--time-limit', type=float, default=60)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology, NAMESPACES
    from taxon_consistency import adjust_ortools as tc
    from analysis.evaluate_lineage_inference import (parse_ncbitaxon_obo, get_true_lineage,
                                                     evaluate_single_lineage, RANKS_7)
    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    hierarchy = defaultdict(set, {g: set(obj['is_a']) for g, obj in go.ont.items()})
    raw = go.normalize_scores(tc.load_predictions(a.predictions), NAMESPACES[a.aspect], {})
    raw = go.propagate_scores(raw, {})
    in_tax, never_tax = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    in_tax, never_tax = go.normalize_constraints(in_tax), go.normalize_constraints(never_tax)
    tax_isa, tax_dis, tax_union = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    parent_map, rank_map, _ = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
    ranked = {t for t, r in rank_map.items() if r in RANKS_7}
    true_lineage = get_true_lineage(a.taxon, parent_map)

    cds = []
    with open(a.cds) as fh:
        rd = csv.reader(fh, delimiter='\t')
        next(rd, None)
        for row in rd:
            cds.append((int(row[0]), int(row[1]), row[3]))
    rng = random.Random(a.seed)
    rows = []
    for length in map(int, a.lengths.split(',')):
        for i in range(a.samples):
            start = rng.randint(1, max(1, a.genome_size - length))
            end = start + length
            prots = sorted({n for s, e, n in cds if s <= end and e >= start and n in raw})
            contig = {p: raw[p] for p in prots}
            ann = {p: {g for g, s in sc.items() if s > a.tau} for p, sc in contig.items()}
            row = {'length': length, 'sample': i, 'start': start, 'n_proteins': len(prots),
                   'n_active': sum(map(len, ann.values()))}
            if not row['n_active']:
                row['category'] = 'no_annotations'
                rows.append(row)
                continue
            costs, _ = tc.compute_demotion_flip_cost(contig, a.tau)
            unc = {}
            t0 = time.monotonic()
            try:
                _, taxa, flips = tc.solve_sat_ortools_taxon_assignment(
                    ann, in_tax, never_tax, hierarchy, costs, tax_isa, tax_dis, tax_union,
                    uncertainty=unc, query_taxa=ranked, time_limit=a.time_limit)
            except RuntimeError as e:
                row.update(category='solver_error', error=str(e)[:200])
                rows.append(row)
                continue
            row['seconds'] = time.monotonic() - t0
            row['status'] = unc.get('status')
            required = set(unc.get('required', []))
            row['n_undetermined'] = len(unc.get('undetermined', []))
            if unc.get('status') != 'OPTIMAL':
                row['category'] = 'not_certified'
            elif not required:
                row['category'] = 'unresolved'
            else:
                ev = evaluate_single_lineage(true_lineage, required, rank_map, parent_map)
                row.update(deepest=ev['deepest_correct_rank'], divergence=ev['first_divergence_rank'])
                row['category'] = 'wrong' if not ev['is_exact_prefix'] else 'correct_' + ev['deepest_correct_rank']
            rows.append(row)
        c = Counter(r['category'] for r in rows if r['length'] == length)
        print(length, dict(c), flush=True)
    summary = defaultdict(Counter)
    for r in rows:
        summary[r['length']][r['category']] += 1
    Path(a.output).write_text(json.dumps({'taxon': a.taxon, 'aspect': a.aspect, 'tau': a.tau,
                                          'predictions': a.predictions, 'genome_size': a.genome_size,
                                          'summary': {str(k): dict(v) for k, v in summary.items()},
                                          'rows': rows}))


if __name__ == '__main__':
    main()
