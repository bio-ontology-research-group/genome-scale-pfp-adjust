"""Evidence-restricted, abstaining complex repair.

Evidence maps use (GO term, protein) keys and nonempty source descriptions.
They are authorizations supplied by the caller, not evidence inferred from a
missing score. Positive evidence propagates to ancestors; negative evidence
propagates to descendants. No demotions are authorized by default.
"""
import math
from collections import Counter, defaultdict
from . import adjust_ortools as cc


def repair(predictions, hierarchy, complex_terms, promotion_evidence,
           demotion_evidence=None, protected_pairs=(), threshold=0.001,
           singleton_penalty=1.0, time_limit=120):
    if not 0 <= threshold < 1 or not math.isfinite(singleton_penalty) or singleton_penalty <= 0:
        raise ValueError('Require threshold in [0,1) and positive finite singleton penalty')
    if round(singleton_penalty * 1000000) < 1:
        raise ValueError('Singleton penalty must be positive at the integer objective scale')
    for protein, scores in predictions.items():
        if any(not math.isfinite(s) or not 0 <= s <= 1 for s in scores.values()):
            raise ValueError(f'Invalid scores for {protein}')
    active, _ = cc.get_annotated_predictions(predictions, threshold)
    if any(parent not in active[p] for p, terms in active.items()
           for g in terms for parent in hierarchy.get(g, ())):
        raise ValueError('Conservative repair requires ancestor-closed active inputs')
    protected = set(protected_pairs)
    if any(p not in active or g not in active[p] for g, p in protected):
        raise ValueError('Protected pairs must be active in the input')
    authorized = []
    for evidence, expand in [(promotion_evidence, cc.get_ancestors),
                             (demotion_evidence or {}, cc.get_descendants)]:
        pairs = set()
        for (g, p), source in evidence.items():
            if p not in predictions or not isinstance(source, str) or not source.strip():
                raise ValueError('Evidence must name an input protein and a nonempty source')
            pairs.update((term, p) for term in {g} | expand(g, hierarchy))
        authorized.append(pairs)
    promotions, demotions = authorized
    protected = {(a, p) for g, p in protected for a in {g} | cc.get_ancestors(g, hierarchy)}
    if demotions & protected:
        raise ValueError('Demotion evidence conflicts with protected annotations')
    # Seed every allowed edit, including non-complex evidence whose ancestors
    # can affect a complex. No top-k restriction is imposed on evidence.
    seeds = defaultdict(set)
    for g, p in promotions | demotions:
        seeds[g].add(p)
    counts = Counter(g for terms in active.values() for g in terms & complex_terms)
    for p, terms in active.items():
        for g in terms & complex_terms:
            if counts[g] == 1:
                seeds[g].add(p)
    pairs, affected = cc.prepare_sparse_domain(active, hierarchy, seeds, complex_terms)
    costs, _ = cc.compute_flip_cost(predictions, threshold)
    stats = {}
    decisions, flips = cc.solve_sat_ortools_optimized(
        active, pairs, hierarchy, costs, affected,
        violation_penalty=singleton_penalty, time_limit=time_limit,
        solver_stats=stats, allowed_promotions=promotions,
        allowed_demotions=demotions, protected_pairs=protected)
    final = cc.get_final_predictions(decisions, predictions, set(), threshold)
    after, _ = cc.get_annotated_predictions(final, threshold)
    added = {(g, p) for p in active for g in after[p] - active[p]}
    removed = {(g, p) for p in active for g in active[p] - after[p]}
    if not added <= promotions or not removed <= demotions or removed & protected:
        raise AssertionError('Export violates evidence restrictions')
    if any(parent not in after[p] for p, terms in after.items()
           for g in terms for parent in hierarchy.get(g, ())):
        raise AssertionError('Export violates GO hierarchy')
    unresolved = []
    for g in sorted(complex_terms):
        members = sorted(p for p, terms in after.items() if g in terms)
        if len(members) == 1:
            candidates = sorted(p for term, p in promotions if term == g and g not in active[p])
            can_remove = (g, members[0]) in demotions and (g, members[0]) not in protected
            reason = ('no_authorized_edit' if not candidates and not can_remove
                      else 'objective_or_coupled_constraints' if stats['optimal']
                      else 'search_incomplete')
            unresolved.append(dict(term=g, protein=members[0], reason=reason,
                                   authorized_candidates=len(candidates)))
    return final, dict(solver=stats, flips=flips,
                      promotions=[list(x) for x in sorted(added)],
                      demotions=[list(x) for x in sorted(removed)],
                      unresolved=unresolved, singleton_penalty=singleton_penalty,
                      threshold=threshold)
