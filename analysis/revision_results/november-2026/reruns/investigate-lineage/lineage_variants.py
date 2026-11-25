"""All-optimum lineage inference under alternative lineage-inference objectives.

Stage 1 without a taxon is a term-level model: one Boolean per GO term, cost c_j
for demoting term j. A variant defines (active terms A, costs c). A ranked taxon
counts as inferred only if every optimal repair requires it (the solver's
feasibility check of its negation at the proven optimum), exactly as in the
revised lineage evaluation (analysis/revision_extensions.py lineage).

Variants (all computed from termstats.py summaries of the canonical,
ancestor-propagated scores):
  margin          A = {j: max_j > tau}, c_j = sum_i (s_ij - tau)          [paper baseline]
  uniform         A as margin, c_j = 1                                    [count each term once]
  max             A as margin, c_j = max_j - tau                          [evidence strength]
  ic              A as margin, c_j = IC_j * sum_i (s_ij - tau)            [IC-weighted mass]
  icmax           A as margin, c_j = IC_j * (max_j - tau)
  thr<t>          A = {j: max_j >= t}, c_j = sum_i (s_ij - t)_+            [high-evidence threshold]
  maxthr<t>       A = {j: max_j >= t}, c_j = max_j
  cnt<t>          A = {j: max_j >= t}, c_j = #proteins with s_ij >= t
  sig<variant>    as <variant>, each constraint signature's terms share its cost (mean over its terms)
  root<variant>   as <variant>, but only terms carrying a taxon constraint cost (others EPS)
  logit[<t>]      A = {j: max_j > t}, c_j = logit(max_j) - logit(t)  [MAP taxon, P(present)=max_j; t=0.5]
"""
import argparse
import gzip
import json
import math
import sys
import time
from collections import defaultdict
from pathlib import Path

EPS = 1e-6


def summarize(st, edges, t):
    """Return (count, margin sum) over scores >= t; t must be a bin edge (or below edge 1)."""
    n = m = 0.0
    for b, e in enumerate(edges):
        if e >= t - 1e-12 or (b == 0 and t <= edges[1]):
            n += st['n'][b]
            m += st['s'][b] - t * st['n'][b]
    return n, m


HIST_EDGES = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9]


def hist(values):
    h = [0] * (len(HIST_EDGES) + 1)
    for v in values:
        h[sum(v >= e for e in HIST_EDGES)] += 1
    return h


def variant_inputs(name, stats, tau, ic, constrained=None):
    if name.startswith('sig'):
        # Each distinct constraint signature (only_in set, never_in set) counts
        # once: its active constrained terms share the signature's evidence
        # (cost / number of active terms with that signature), so a clade's
        # evidence does not grow with the size of its GO vocabulary.
        act, cost = variant_inputs(name[3:], stats, tau, ic)
        groups = defaultdict(list)
        for g in act:
            if g in constrained:
                groups[constrained[g]].append(g)
        out = {g: EPS for g in act}
        for members in groups.values():
            for g in members:
                out[g] = max(cost[g] / len(members), EPS)
        return act, out
    if name.startswith('root'):
        # Count evidence once per taxon constraint: only terms that carry an
        # only_in/never_in constraint cost anything; their unconstrained
        # descendants (already summarized by max-propagation) are free.
        act, cost = variant_inputs(name[4:], stats, tau, ic)
        return act, {g: (c if g in constrained else EPS) for g, c in cost.items()}
    edges = stats['bin_edges']
    terms = stats['terms']
    act, cost = set(), {}
    for g, st in terms.items():
        mx = st['max']
        if name.startswith('logit'):
            # MAP taxon under independent term presence with P(j present) = max_j:
            # keeping j gains log p, demoting gains log(1-p); terms with p <= 0.5
            # are better absent (free to demote), so A = {p > 0.5}, c = logit(p).
            t = float(name[5:]) if len(name) > 5 else 0.5
            if mx <= t:
                continue
            p = min(mx, 1 - 1e-3)
            c = math.log(p / (1 - p)) - math.log(t / (1 - t))
            act.add(g)
            cost[g] = max(c, EPS)
            continue
        if name in ('margin', 'uniform', 'max', 'ic', 'icmax'):
            if mx <= tau:
                continue
            m = stats['margin'][g]
            c = {'margin': m, 'uniform': 1.0, 'max': mx - tau,
                 'ic': ic.get(g, 1.0) * m, 'icmax': ic.get(g, 1.0) * (mx - tau)}[name]
        else:
            kind, t = (name[:6], float(name[6:])) if name.startswith('maxthr') else (name[:3], float(name[3:]))
            if mx < t:
                continue
            n, m = summarize(st, edges, t)
            c = {'thr': m, 'maxthr': mx, 'cnt': n}[kind]
        act.add(g)
        cost[g] = max(c, EPS)
    return act, cost


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--resources', required=True)
    ap.add_argument('--termstats', required=True, help='dir with <pred>/<aspect>/<taxon>.json.gz')
    ap.add_argument('--predictor', required=True)
    ap.add_argument('--taxa', required=True, help='comma-separated taxa')
    ap.add_argument('--variants', required=True)
    ap.add_argument('--aspects', default='cc,mf,bp,joint')
    ap.add_argument('--ic-file')
    ap.add_argument('--truth', action='store_true', help='also solve with the true taxon fixed')
    ap.add_argument('--truth-variants', default='margin,thr0.5,logit,sigthr0.5',
                    help='variants for which the true-taxon solve is run (with --truth)')
    ap.add_argument('--time-limit', type=float, default=120)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from gaef_patches.ontology import Ontology
    from taxon_consistency import adjust_ortools as tc
    from analysis.evaluate_lineage_inference import (parse_ncbitaxon_obo, get_true_lineage,
                                                     evaluate_single_lineage, RANKS_7)
    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    hierarchy = {g: set(o['is_a']) for g, o in go.ont.items()}
    only, never = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    only, never = go.normalize_constraints(only), go.normalize_constraints(never)
    isa, dis, uni = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    parents, ranks, names = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
    query = {t for t, r in ranks.items() if r in RANKS_7}
    constrained = {g: (frozenset(only.get(g, ())), frozenset(never.get(g, ()))) for g in set(only) | set(never)}
    ic = {}
    if a.ic_file:
        for line in open(a.ic_file):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2:
                try:
                    ic[go.canonical_id(p[0]) or p[0]] = float(p[1])
                except ValueError:
                    pass
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for taxon in a.taxa.split(','):
        target = out / f'{a.predictor}_{taxon}.json'
        done = json.loads(target.read_text()) if target.exists() else {'records': []}
        have = {(r['variant'], r['aspect']) for r in done['records']}
        stats = {}
        for asp in ['cc', 'mf', 'bp']:
            with gzip.open(Path(a.termstats) / a.predictor / asp / f'{taxon}.json.gz', 'rt') as fh:
                stats[asp] = json.load(fh)
        truth = get_true_lineage(taxon, parents)
        for var in a.variants.split(','):
            per = {asp: variant_inputs(var, stats[asp], stats[asp]['tau'], ic, constrained) for asp in stats}
            for asp in a.aspects.split(','):
                if (var, asp) in have:
                    continue
                if asp == 'joint':
                    act = set().union(*(per[x][0] for x in per))
                    cost = {}
                    for x in per:
                        for g, c in per[x][1].items():
                            cost[g] = cost.get(g, 0.0) + c
                else:
                    act, cost = per[asp]
                rec = {'variant': var, 'aspect': asp, 'taxon': taxon, 'n_active': len(act)}
                if not act:
                    rec.update(status='EMPTY', required=[], category='unresolved')
                    done['records'].append(rec)
                    continue
                unc = {}
                t0 = time.monotonic()
                try:
                    adj, assign, _ = tc.solve_sat_ortools_taxon_assignment(
                        {'genome': act}, only, never, hierarchy, cost, isa, dis, uni,
                        uncertainty=unc, query_taxa=query, time_limit=a.time_limit)
                except RuntimeError as e:
                    rec.update(status='ERROR', error=str(e)[:200], required=[], category='not_certified')
                    done['records'].append(rec)
                    continue
                rec['seconds'] = round(time.monotonic() - t0, 2)
                rec['status'] = unc.get('status')
                req = sorted(unc.get('required', []))
                rec['required'] = req
                rec['possible'] = sorted(unc.get('possible', []))
                rec['n_undetermined'] = len(unc.get('undetermined', []))
                demoted = sorted(g for g, keep in adj['genome'].items() if not keep)
                rec['n_demoted'] = len(demoted)
                rec['demoted_cost'] = sum(cost.get(g, 0) for g in demoted)
                top = sorted(demoted, key=lambda g: -cost.get(g, 0))[:15]
                rec['top_demoted'] = [(g, round(cost.get(g, 0), 4)) for g in top]
                mx = {g: st['max'] for x in stats for g, st in stats[x]['terms'].items()}
                rec['demoted_max_hist'] = hist([mx.get(g, 0) for g in demoted if g in act])
                if a.truth and var in a.truth_variants.split(','):
                    try:
                        tadj, _, _ = tc.solve_sat_ortools_taxon_assignment(
                            {'genome': act}, only, never, hierarchy, cost, isa, dis, uni,
                            genome_taxon_id=taxon, time_limit=a.time_limit)
                        tdem = sorted(g for g, keep in tadj['genome'].items() if not keep)
                        rec['truth_n_demoted'] = len(tdem)
                        rec['truth_cost'] = sum(cost.get(g, 0) for g in tdem)
                        rec['truth_top_demoted'] = [(g, round(cost.get(g, 0), 4), round(mx.get(g, 0), 3))
                                                    for g in sorted(tdem, key=lambda g: -cost.get(g, 0))[:15]]
                        rec['truth_demoted_max_hist'] = hist([mx.get(g, 0) for g in tdem if g in act])
                        rec['truth_only_demoted'] = len(set(tdem) - set(demoted))
                        rec['opt_only_demoted'] = len(set(demoted) - set(tdem))
                    except RuntimeError as e:
                        rec['truth_error'] = str(e)[:200]
                ev = evaluate_single_lineage(truth, set(req), ranks, parents)
                rec['deepest'] = ev['deepest_correct_rank']
                rec['divergence'] = ev['first_divergence_rank']
                rec['exact_prefix'] = ev['is_exact_prefix']
                rec['inferred_ranks'] = ev['inferred_ranks']
                if rec['status'] != 'OPTIMAL':
                    rec['category'] = 'not_certified'
                elif not ev['inferred_ranks']:
                    rec['category'] = 'unresolved'
                elif not ev['is_exact_prefix']:
                    rec['category'] = 'wrong'
                else:
                    rec['category'] = 'correct'
                done['records'].append(rec)
                print(a.predictor, taxon, var, asp, rec['category'], rec['deepest'], rec['divergence'],
                      rec.get('seconds'), flush=True)
            target.write_text(json.dumps(done))
        target.write_text(json.dumps(done))


if __name__ == '__main__':
    main()
