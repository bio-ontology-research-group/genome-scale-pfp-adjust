"""Contig lineage experiment (as todo-reruns-2026-10-31/contig_lineage.py) under a
lineage-inference objective variant from lineage_variants.py.

Same sampling (seed 42, contig lengths in order), same proteins, same all-optimum
required-taxon evaluation; only the term set and costs fed to Stage 1 change.
"""
import argparse
import csv
import json
import random
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from termstats import BIN_EDGES  # noqa: E402
from lineage_variants import variant_inputs  # noqa: E402


def contig_stats(contig, tau):
    nb = len(BIN_EDGES)
    terms = {}
    for sc in contig.values():
        for g, s in sc.items():
            if s <= tau:
                continue
            st = terms.get(g)
            if st is None:
                st = terms[g] = {'n': [0] * nb, 's': [0.0] * nb, 'max': 0.0}
            b = 0
            while b + 1 < nb and s >= BIN_EDGES[b + 1]:
                b += 1
            st['n'][b] += 1
            st['s'][b] += s
            st['max'] = max(st['max'], s)
    margin = {g: sum(st['s']) - tau * sum(st['n']) for g, st in terms.items()}
    return {'tau': tau, 'bin_edges': BIN_EDGES, 'terms': terms, 'margin': margin}


def main():
    ap = argparse.ArgumentParser()
    for k in ['code', 'resources', 'predictions', 'cds', 'taxon', 'aspect', 'variants', 'output']:
        ap.add_argument('--' + k, required=True)
    ap.add_argument('--genome-size', type=int, required=True)
    ap.add_argument('--tau', type=float, required=True)
    ap.add_argument('--lengths', default='1000,10000,100000,1000000')
    ap.add_argument('--keep-lengths', default='100000,1000000')
    ap.add_argument('--samples', type=int, default=1000)
    ap.add_argument('--evaluate-first', type=int, default=1000,
                    help='evaluate only the first K samples per length (RNG stream unchanged, so these are '
                         'the same contigs as samples 0..K-1 of the original 1000-sample run)')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--time-limit', type=float, default=60)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology, NAMESPACES
    from taxon_consistency import adjust_ortools as tc
    from analysis.evaluate_lineage_inference import (parse_ncbitaxon_obo, get_true_lineage,
                                                     evaluate_single_lineage, RANKS_7)
    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    hierarchy = defaultdict(set, {g: set(obj['is_a']) for g, obj in go.ont.items()})
    raw = go.propagate_scores(go.normalize_scores(tc.load_predictions(a.predictions), NAMESPACES[a.aspect], {}), {})
    in_tax, never_tax = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    in_tax, never_tax = go.normalize_constraints(in_tax), go.normalize_constraints(never_tax)
    tax_isa, tax_dis, tax_union = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    parent_map, rank_map, _ = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
    ranked = {t for t, r in rank_map.items() if r in RANKS_7}
    true_lineage = get_true_lineage(a.taxon, parent_map)
    constrained = {g: (frozenset(in_tax.get(g, ())), frozenset(never_tax.get(g, ()))) for g in set(in_tax) | set(never_tax)}
    cds = []
    with open(a.cds) as fh:
        rd = csv.reader(fh, delimiter='\t')
        next(rd, None)
        for row in rd:
            cds.append((int(row[0]), int(row[1]), row[3]))
    keep = set(map(int, a.keep_lengths.split(',')))
    variants = a.variants.split(',')
    rng = random.Random(a.seed)
    rows = []
    for length in map(int, a.lengths.split(',')):
        for i in range(a.samples):
            start = rng.randint(1, max(1, a.genome_size - length))  # identical RNG stream
            if length not in keep or i >= a.evaluate_first:
                continue
            end = start + length
            prots = sorted({n for s, e, n in cds if s <= end and e >= start and n in raw})
            stats = contig_stats({p: raw[p] for p in prots}, a.tau)
            for var in variants:
                act, cost = variant_inputs(var, stats, a.tau, {}, constrained)
                row = {'variant': var, 'length': length, 'sample': i, 'start': start,
                       'n_proteins': len(prots), 'n_active': len(act)}
                if not act:
                    row['category'] = 'no_annotations'
                    rows.append(row)
                    continue
                unc = {}
                t0 = time.monotonic()
                try:
                    tc.solve_sat_ortools_taxon_assignment(
                        {'contig': act}, in_tax, never_tax, hierarchy, cost, tax_isa, tax_dis, tax_union,
                        uncertainty=unc, query_taxa=ranked, time_limit=a.time_limit)
                except RuntimeError as e:
                    row.update(category='solver_error', error=str(e)[:200])
                    rows.append(row)
                    continue
                row['seconds'] = round(time.monotonic() - t0, 2)
                required = set(unc.get('required', []))
                ev = evaluate_single_lineage(true_lineage, required, rank_map, parent_map)
                row['required'] = sorted(required)
                if unc.get('status') != 'OPTIMAL':
                    row['category'] = 'not_certified'
                elif not ev['inferred_ranks']:
                    row['category'] = 'unresolved'
                else:
                    row.update(deepest=ev['deepest_correct_rank'], divergence=ev['first_divergence_rank'])
                    row['category'] = 'wrong' if not ev['is_exact_prefix'] else 'correct_' + ev['deepest_correct_rank']
                rows.append(row)
        for var in variants:
            c = Counter(r['category'] for r in rows if r['length'] == length and r['variant'] == var)
            if c:
                print(length, var, dict(c), flush=True)
    summary = defaultdict(lambda: defaultdict(Counter))
    for r in rows:
        summary[r['variant']][r['length']][r['category']] += 1
    Path(a.output).write_text(json.dumps({
        'taxon': a.taxon, 'aspect': a.aspect, 'tau': a.tau, 'predictions': a.predictions,
        'summary': {v: {str(k): dict(c) for k, c in d.items()} for v, d in summary.items()}, 'rows': rows}))


if __name__ == '__main__':
    main()
