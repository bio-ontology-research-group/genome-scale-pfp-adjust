"""Alternative Stage 2 policies on one predictor x organism (Stage 1 input).

All policies keep H1 (singletons only) and H2 (sparse closed domain) and use
the existing solver; they differ only in which edits are authorized and, for
one variant, in the flip cost. Every policy uses the soft singleton
formulation with a large penalty (P=100 >> any repair cost), so a singleton is
repaired whenever an authorized edit exists and is otherwise left (abstain).

Label-free promotion evidence for (X, q), never derived from curated sets:
  score : primary score s_q,X > floor (floor=0; DeepGO-SE: s >= --dense-floor)
  rel   : q active (Stage 1) for a direct is_a parent of X that is itself an
          operational complex term other than GO:0032991
  xpred : at least one other predictor's Stage 1 output has X active on q
Policies:
  top50        H3 candidates (top-50 by score), promotion or demotion (= standard, soft)
  top50_nodem  H3 candidates, no demotion
  logodds      H3 candidates, cost |logit(s)-logit(tau)| (s clipped to [1e-6,1-1e-6])
  ev_hard      evidence candidates (<=50 by evidence count, score), demotion allowed
  ev_prefer    evidence candidates; demotion only for singletons without evidence
  ev_abstain   evidence candidates; no demotion (abstain without evidence)
"""
import argparse, gzip, json, math, time
from collections import defaultdict
from pathlib import Path
import common_ic as C

GENERIC = {'GO:0032991'}


def load_selected(path, terms):
    out = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            d = {}
            for tok in cols[1:]:
                g, s = tok.split('|')
                if g in terms:
                    d[g] = float(s)
            out[cols[0]] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ['code', 'go', 'accepted', 'predictor', 'taxon', 'out']:
        ap.add_argument('--'+k, required=True)
    ap.add_argument('--variants', default='acc_standard,acc_demotion,top50_nodem,logodds,ev_hard,ev_prefer,ev_abstain,sx_hard,sx_prefer,sx_abstain')
    ap.add_argument('--dense-floor', type=float, default=None,
                    help='DeepGO-SE score floor for score evidence (None: no score evidence for DeepGO-SE)')
    ap.add_argument('--time-limit', type=float, default=120)
    ap.add_argument('--penalty', type=float, default=100.)
    ap.add_argument('--workers', type=int, default=1)
    a = ap.parse_args()
    go, H, isa_h, complexes, cc = C.resources(a.code, a.go)
    if a.workers > 1:
        class Solver(cc.cp_model.CpSolver):
            def Solve(self, model, *x, **k):
                self.parameters.num_search_workers = a.workers
                return super().Solve(model, *x, **k)
        cc.cp_model.CpSolver = Solver
    tau = C.TAU[a.predictor]
    s1 = C.stage_file(a.accepted, a.predictor, 'stage1', a.taxon)
    scores = cc.load_predictions(s1)
    ann, terms = cc.get_annotated_predictions(scores, tau)
    members = defaultdict(set)
    for p, ts in ann.items():
        for g in ts & complexes:
            members[g].add(p)
    singles = sorted(g for g, m in members.items() if len(m) == 1)
    lone = {X: next(iter(members[X])) for X in singles}
    names = sorted(scores)
    anc_c, desc_c = {}, {}
    def anc(g):
        if g not in anc_c: anc_c[g] = cc.get_ancestors(g, H)
        return anc_c[g]
    def desc(g):
        if g not in desc_c: desc_c[g] = cc.get_descendants(g, H)
        return desc_c[g]
    # --- evidence ---
    others = {}
    for op in C.PREDS:
        if op == a.predictor:
            continue
        sc = load_selected(C.stage_file(a.accepted, op, 'stage1', a.taxon), set(singles))
        others[op] = {(g, p) for p, d in sc.items() for g, v in d.items() if v > C.TAU[op]}
        del sc
    ev = {}
    top50 = {}
    for X in singles:
        L = lone[X]
        order = sorted(names, key=lambda p: (-scores[p].get(X, 0.0), p))
        top50[X] = [p for p in order[:50] if p != L]
        parents = {y for y in isa_h.get(X, ()) if y in complexes and y not in GENERIC}
        e = defaultdict(list)
        for p in names:
            if p == L:
                continue
            s = scores[p].get(X, 0.0)
            if (a.predictor != 'deepgo-se' and s > 0) or (a.predictor == 'deepgo-se' and a.dense_floor is not None and s >= a.dense_floor):
                e[p].append('score')
            if parents & ann[p]:
                e[p].append('rel')
            if any((X, p) in o for o in others.values()):
                e[p].append('xpred')
        ev[X] = {p: t for p, t in e.items()}
    evcand = {X: sorted(ev[X], key=lambda p: (-len(ev[X][p]), -scores[p].get(X, 0.0), p))[:50] for X in singles}
    # strict evidence: primary score or other-predictor support (no parent-complex 'rel' evidence)
    sxcand = {X: sorted((p for p, t in ev[X].items() if set(t) & {'score', 'xpred'}),
                        key=lambda p: (-len(ev[X][p]), -scores[p].get(X, 0.0), p))[:50] for X in singles}
    base_costs, _ = cc.compute_flip_cost(scores, tau)

    def logodds_costs():
        lt = math.log(tau/(1-tau))
        f = defaultdict(lambda: defaultdict(lambda: abs(math.log(1e-6/(1-1e-6)) - lt)))
        for p, d in scores.items():
            for g, s in d.items():
                s = min(max(s, 1e-6), 1-1e-6)
                f[p][g] = abs(math.log(s/(1-s)) - lt)
        return f

    def run(cand, demotable, costs, hard=False):
        part = {X: sorted({lone[X]} | set(cand[X])) for X in singles}
        pairs, affected = cc.prepare_sparse_domain(ann, H, part, complexes)
        prom = {(g, q) for X in singles for q in cand[X] for g in {X} | anc(X) if g not in ann[q]}
        dem = {(g, lone[X]) for X in demotable for g in {X} | desc(X) if g in ann[lone[X]]}
        st = {}
        t0 = time.monotonic()
        try:
            sol, flips = cc.solve_sat_ortools_optimized(ann, pairs, H, costs, affected, violation_penalty=None if hard else a.penalty,
                                                        time_limit=a.time_limit, solver_stats=st,
                                                        allowed_promotions=prom, allowed_demotions=dem)
            st['formulation'] = 'hard' if hard else 'soft'
        except RuntimeError:
            st = {}
            sol, flips = cc.solve_sat_ortools_optimized(ann, pairs, H, costs, affected, violation_penalty=a.penalty,
                                                        time_limit=a.time_limit, solver_stats=st,
                                                        allowed_promotions=prom, allowed_demotions=dem)
            st['formulation'] = 'soft_fallback'
        st['wall'] = time.monotonic() - t0
        return cc.get_final_predictions(sol, scores, terms, tau), st

    out = dict(predictor=a.predictor, taxon=a.taxon, tau=tau, n_singletons=len(singles), dense_floor=a.dense_floor,
               evidence={X: ev[X] for X in singles}, lone=lone,
               lone_score={X: scores[lone[X]].get(X, 0.0) for X in singles}, variants={})
    nodem = set()
    for v in a.variants.split(','):
        if v in ('acc_standard', 'acc_demotion'):
            final = cc.load_predictions(C.stage_file(a.accepted, a.predictor, 'stage2_' + v[4:], a.taxon))
            st = {'source': 'accepted'}
        elif v == 'top50':
            final, st = run(top50, set(singles), base_costs)
        elif v == 'top50_nodem':
            final, st = run(top50, nodem, base_costs, hard=True)
        elif v == 'logodds':
            final, st = run(top50, set(singles), logodds_costs(), hard=True)
        elif v == 'ev_hard':
            final, st = run(evcand, set(singles), base_costs, hard=True)
        elif v == 'ev_prefer':
            final, st = run(evcand, {X for X in singles if not evcand[X]}, base_costs, hard=True)
        elif v == 'sx_hard':
            final, st = run(sxcand, set(singles), base_costs, hard=True)
        elif v == 'sx_prefer':
            final, st = run(sxcand, {X for X in singles if not sxcand[X]}, base_costs, hard=True)
        elif v == 'sx_abstain':
            final, st = run(sxcand, nodem, base_costs)
        elif v == 'ev_abstain':
            final, st = run(evcand, nodem, base_costs)
        else:
            raise ValueError(v)
        after = C.active(final, tau)
        prom = sorted([g, p] for p in after for g in after[p] - ann.get(p, set()))
        dem = sorted([g, p] for p in ann for g in ann[p] - after.get(p, set()))
        cnt = defaultdict(int)
        for p, ts in after.items():
            for g in ts & complexes:
                cnt[g] += 1
        present = [g for g in cnt if cnt[g]]
        sing_after = [g for g in present if cnt[g] == 1]
        outcome = {}
        for X in singles:
            L = lone[X]
            if X not in after.get(L, set()):
                outcome[X] = 'demoted'
            elif cnt[X] >= 2:
                outcome[X] = 'promoted'
            else:
                outcome[X] = 'unresolved'
        if not v.startswith('acc_'):
            dest = Path(a.out)/a.predictor/v
            dest.mkdir(parents=True, exist_ok=True)
            cc.save_predictions(final, dest/s1.name)
        out['variants'][v] = dict(solver=st, promotions=prom, demotions=dem, outcome=outcome,
                                  n_present=len(present), n_singletons_after=len(sing_after))
        print(v, st, len(prom), len(dem), len(sing_after), flush=True)
        del final, after
    with gzip.open(Path(a.out)/a.predictor/f'{a.taxon}.json.gz', 'wt') as fh:
        json.dump(out, fh)


if __name__ == '__main__':
    main()
