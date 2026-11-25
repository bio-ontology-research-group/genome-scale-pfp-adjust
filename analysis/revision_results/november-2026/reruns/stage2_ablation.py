"""Stage 2 reduction/heuristic ablation with the corrected implementation.

Input: a known-taxon Stage 1 output (canonical, ancestor-propagated) of one
organism. Configurations:
  naive      every CC term x every protein, cardinality constraint on every complex term
  H1+H2      incoherent complexes only, sparse domain closed over affected constraints,
             every protein a promotion candidate
  H1+H2+H3   as H1+H2 with the 50 highest-scoring promotion candidates per complex
Naive variable/constraint counts are computed analytically when the model is
not built. Objectives are compared to check that the reductions preserve the
optimum where the naive model can be solved.
"""
import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--code', required=True)
    ap.add_argument('--go', required=True)
    ap.add_argument('--stage1', required=True)
    ap.add_argument('--taxon', required=True)
    ap.add_argument('--tau', type=float, default=.001)
    ap.add_argument('--timeout', type=float, default=600)
    ap.add_argument('--solve-naive', action='store_true')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--output', required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    from ortools.sat.python import cp_model
    from gaef_patches.ontology import Ontology
    from complex_coherence import adjust_ortools as cc

    captured = {}

    class Model(cp_model.CpModel):
        def __init__(self, *x, **k):
            super().__init__(*x, **k)
            captured['model'] = self

    class Solver(cp_model.CpSolver):
        def Solve(self, model, *x, **k):
            self.parameters.max_time_in_seconds = a.timeout
            self.parameters.num_search_workers = a.workers
            t = time.monotonic()
            status = super().Solve(model, *x, **k)
            captured['solve_s'] = time.monotonic() - t
            captured['status'] = self.StatusName(status)
            captured['objective'] = self.ObjectiveValue() / 1e6 if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
            captured['bound'] = self.BestObjectiveBound() / 1e6
            return status

    cc.cp_model.CpModel = Model
    cc.cp_model.CpSolver = Solver

    go = Ontology(a.go)
    cc_terms = go.get_namespace_terms('cellular_component')
    hierarchy = defaultdict(set, {g: {p for p in obj['is_a'] if p in cc_terms}
                                  for g, obj in go.ont.items() if g in cc_terms})
    isa_go = Ontology(a.go, with_rels=False)
    isa_h = {g: set(obj['is_a']) for g, obj in isa_go.ont.items()}
    homomers = go.normalize_terms(cc.load_homodimer_terms(Path(a.code) / 'complex_coherence/protein_complexes.tsv'))
    complexes = cc.get_heteromeric_complexes(homomers, set(), isa_h)
    raw = cc.load_predictions(a.stage1)
    raw = {p: {g: s for g, s in sc.items() if g in cc_terms} for p, sc in raw.items()}
    ann, terms = cc.get_annotated_predictions(raw, a.tau)
    costs, _ = cc.compute_flip_cost(raw, a.tau)
    n_p = len(ann)
    counts = defaultdict(int)
    for ts in ann.values():
        for g in ts & complexes:
            counts[g] += 1
    incoherent = sorted(g for g, n in counts.items() if n == 1)
    out = {'taxon': a.taxon, 'tau': a.tau, 'n_proteins': n_p, 'n_active_pairs': sum(map(len, ann.values())),
           'n_incoherent_complexes': len(incoherent), 'configs': {}}
    # naive: no H1, H2 or H3. Every CC term x every protein is a decision and every
    # complex term carries the cardinality constraint (same encoding as the sparse solver).
    all_terms = sorted(cc_terms)
    all_complexes = sorted(complexes & cc_terms)
    edges = sum(len(hierarchy.get(g, ())) for g in all_terms)
    out['configs']['naive'] = {'n_vars': len(all_terms) * n_p + 2 * len(all_complexes),
                               'n_constraints': edges * n_p + 2 * len(all_complexes),
                               'analytic': True}
    if a.solve_naive:
        captured.clear()
        t = time.monotonic()
        pairs = {(g, p) for g in all_terms for p in ann}
        cand = {g: sorted(ann) for g in all_complexes}
        stats = {}
        try:
            sol, flips = cc.solve_sat_ortools_optimized(ann, pairs, hierarchy, costs, cand,
                                                        time_limit=a.timeout, solver_stats=stats)
        except RuntimeError as e:
            flips = None
        proto = captured['model'].Proto()
        out['configs']['naive'].update(n_vars=len(proto.variables), n_constraints=len(proto.constraints),
                                       analytic=False, total_s=time.monotonic() - t, flips=flips,
                                       **{k: captured.get(k) for k in ['solve_s', 'status', 'objective', 'bound']})
        del pairs
        print('naive', out['configs']['naive'], flush=True)
    for label, top_k in [('H1+H2', None), ('H1+H2+H3', 50)]:
        captured.clear()
        t = time.monotonic()
        cand = cc.get_participating_proteins(ann, complexes, raw, top_k=top_k)
        pairs, cand = cc.prepare_sparse_domain(ann, hierarchy, cand, complexes)
        stats = {}
        try:
            sol, flips = cc.solve_sat_ortools_optimized(ann, pairs, hierarchy, costs, cand,
                                                        time_limit=a.timeout, solver_stats=stats)
        except RuntimeError as e:
            flips = None
            stats['status'] = str(e)
        proto = captured['model'].Proto()
        out['configs'][label] = {'n_vars': len(proto.variables), 'n_constraints': len(proto.constraints),
                                 'total_s': time.monotonic() - t, 'flips': flips,
                                 **{k: captured.get(k) for k in ['solve_s', 'status', 'objective', 'bound']}}
        print(label, out['configs'][label], flush=True)
    Path(a.output).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
