"""Rerun one predictor/aspect/organism with explicit inputs and test membership.

COPY of analysis/revision_audit_run.py (commit c260eef) with one addition:
--stage2-gated adds an evidence-gated Stage 2 run ('stage2_gated'), the
'sx_prefer' policy of todo-reruns-2026-10-31/investigate-complex/variants.py:
  * only singleton complexes are repaired (H1), sparse closed domain (H2);
  * protein q != lone member L is a promotion candidate for singleton X only if
    one of the OTHER two timeset predictors has X active on q (score > its native
    tau) in its Stage 1 output (canonicalized, propagated, taxon-repaired);
    candidates ordered by (-#evidence types incl. parent-complex 'rel', -score, id),
    at most 50 per singleton; promotions restricted to X and its ancestors on them;
  * the lone annotation of X is demotable (X and its active descendants on L)
    only if X has no candidate (prefer promotion);
  * hard formulation first; if infeasible, soft fallback with penalty 100.
  (variants.py also admits primary-score evidence s>0 for non-DeepGO predictors;
  it is empty for all 50x3 timeset cases, so the definitions coincide.)
The code manifest uses the pipeline code under $PFP_CODE (ROOT).

Each output is independent and carries input checksums. Aggregation is a
separate step, so missing organisms cannot silently produce a complete table.
"""
import argparse
import os
from collections import defaultdict, Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time
import platform
import os
import importlib.metadata

ROOT = Path(os.environ.get('PFP_CODE', Path(__file__).resolve().parents[1])).resolve()
sys.path.insert(0, str(ROOT))
from analysis.evaluation_scope import load_temporal_proteins, select_annotations, scope_manifest
from gaef_patches.ontology import Ontology, NAMESPACES
from pipeline.genome_scores_evaluator import load_annotations, check_consistency_fast
from taxon_consistency import adjust_ortools as tc
from complex_coherence import adjust_ortools as cc
from gaef_patches import coherence


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def active(raw, tau):
    return {p: {g for g, s in scores.items() if s > tau} for p, scores in raw.items()}


def stats(raw, tau, complexes):
    ann = active(raw, tau)
    counts = Counter(g for terms in ann.values() for g in terms)
    present = set(counts) & complexes
    singletons = sorted(g for g in present if counts[g] == 1)
    return {'n_proteins': len(raw), 'active_predictions': sum(counts.values()),
            'n_active_terms': len(counts),
            'grid_decisions': len(raw) * len(counts),
            'n_present_complexes': len(present), 'n_singleton_complexes': len(singletons),
            'complex_coherence': 1 - len(singletons)/len(present) if present else 1.0}


def run(args):
    if os.environ.get('PFP_COST_MODE', 'margin') != 'margin':
        raise ValueError('This experiment requires margin costs')
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    data, resources = Path(args.data), Path(args.resources)
    files = sorted((data / f'{args.predictor}_results' / args.aspect / 'predictions').glob(f'predictions_fold_*_taxon_{args.taxon}.tsv'))
    if len(files) != 1:
        raise ValueError(f'Expected one prediction file, found {files}')
    pred_file = files[0]
    ann_file = Path(args.annotations or data/'annotations-go-basic') / f'annots_taxon_{args.taxon}.tsv'
    temporal = load_temporal_proteins(args.proteins) if args.proteins else None
    annotations = select_annotations(load_annotations(ann_file), args.taxon, temporal)
    go_file = resources / 'go-basic.obo'
    go = Ontology(go_file)
    hierarchy = defaultdict(set, {g: set(obj['is_a']) for g, obj in go.ont.items()})
    namespace = go.get_namespace_terms(NAMESPACES[args.aspect])
    identifier_audit = {}
    raw = go.normalize_scores(tc.load_predictions(pred_file), NAMESPACES[args.aspect], identifier_audit)
    raw = go.propagate_scores(raw, identifier_audit)
    annotations = {p: go.normalize_terms(ts) for p, ts in annotations.items()}
    constraints_file = resources / 'go_taxon_constraints_extracted_obo.tsv'
    in_tax, never_tax = tc.load_constraints(constraints_file)
    in_tax, never_tax = go.normalize_constraints(in_tax), go.normalize_constraints(never_tax)
    tax_files = [resources / 'taxon_hierarchy.tsv', resources / 'ncbitaxon_hierarchy.tsv']
    tax_isa, tax_disjoint, tax_union = tc.load_taxon_hierarchy(*tax_files)
    lineage = {tc.normalize_taxon_id(args.taxon)}
    lineage.update(tc.get_all_ancestors_of_list(tax_isa, list(lineage)))
    for union, members in tax_union.items():
        if members & lineage:
            lineage.add(union)
    homomers = go.normalize_terms(cc.load_homodimer_terms(ROOT / 'complex_coherence/protein_complexes.tsv'))
    isa_go = Ontology(go_file, with_rels=False)
    isa_hierarchy = {g: set(obj['is_a']) for g, obj in isa_go.ont.items()}
    complexes = cc.get_heteromeric_complexes(homomers, set(), isa_hierarchy) if args.aspect == 'cc' else set()
    cdir = ROOT / 'data/constraints'
    has_part = coherence.parse_has_part(cdir / 'has_part_relations.txt')
    ec2go = coherence.parse_ec2go(cdir / 'ec2go_v2025-03-16')
    pathways = coherence.map_pathways_to_go_terms(cdir / 'metacyc_GO_v2025-03-16_with_EC.tsv', ec2go)
    has_part, ec2go, pathways = coherence.normalize_go_resources(go, has_part, ec2go, pathways)
    ancestor_cache = {}
    native = .3 if args.predictor == 'deepgo-se' else .001
    out = {'schema_version': 2,
           'identifier_normalization': identifier_audit,
           'configuration': {'sweep': args.sweep, 'stage2': args.stage2, 'soft': args.soft, 'infer_taxons': args.infer_taxons, 'solver_seconds': args.solver_seconds, 'top_k': 50,
                             'stage2_gated': bool(args.stage2_gated), 'gated_policy': 'sx_prefer_xpred_stage1_v1' if args.stage2_gated else None,
                             'input_alignment': 'primary_ids_max_score_is_a_part_of_closure_v1'},
           'software': {'python': platform.python_version(), **{pkg: importlib.metadata.version(pkg) for pkg in ['ortools', 'numpy', 'pandas', 'scikit-learn']}}, 'predictor': args.predictor, 'aspect': args.aspect,
           'taxon': args.taxon, 'native_tau': native, 'n_complex_terms': len(complexes),
           'protein_scope': scope_manifest(args.proteins) if args.proteins else {'definition': 'held-out organism'},
           'inputs': {str(f): sha(f) for f in [pred_file, ann_file, go_file, constraints_file, ROOT / 'complex_coherence/protein_complexes.tsv'] + tax_files + sorted(cdir.glob('*')) if f.is_file()},
           'gated_driver_sha256': sha(__file__),
           'code': {str(f.relative_to(ROOT)): sha(f) for folder in ['analysis', 'taxon_consistency', 'complex_coherence', 'gaef_patches', 'pipeline'] for f in sorted((ROOT / folder).rglob('*.py'))},
           'annotations': {p: sorted(ts) for p, ts in annotations.items()},
           'n_mask_proteins': len(temporal.get(args.taxon, ())) if temporal is not None else len(annotations),
           'n_annotation_matches': len(annotations),
           'n_prediction_matches': len(set(annotations) & set(raw)), 'runs': []}
    def record(label, scores, tau, before=None):
        result = {'mode': label, 'tau': tau, **stats(scores, tau, complexes),
                  'test_scores': {p: scores.get(p, {}) for p in annotations}}
        binary = active(scores, tau)
        union = set().union(*binary.values())
        closed = set()
        for term in union:
            if term not in ancestor_cache:
                ancestor_cache[term] = go.get_ancestors(term)
            closed.update(ancestor_cache[term])
        violations = sorted(g for g in closed if
                            any(t not in lineage for t in in_tax.get(g, ())) or
                            any(t in lineage for t in never_tax.get(g, ())))
        result['known_taxon_consistent'] = not violations
        result['gaef_taxon_consistent'] = bool(check_consistency_fast({'genome': closed}, ROOT / 'data/constraints/taxon_constraints.tsv', ontology=go))
        result['violating_terms'] = violations
        process, _ = coherence.check_has_part({'genome': closed}, has_part)
        process_denominator = sum(len(has_part.get(g, ())) for g in closed)
        result['process_coherence'] = process / 100 if process_denominator else None
        result['process_relation_count'] = process_denominator
        if args.aspect in {'bp', 'mf'}:
            _, complete, annotated, _ = coherence.analyze_genome({'genome': closed}, pathways, ec2go)
            result['pathway_coherence'] = len(complete)/len(annotated) if annotated else None
            result['pathway_present'] = len(annotated)
        if args.aspect == 'bp':
            from GAEF.completeness import read_essential_terms
            core, _ = read_essential_terms(cdir / 'essential_terms.tsv')
            core_terms = {go.canonical_id(e['term']) or e['term'] for e in core}
            result['essential_completeness'] = len(core_terms & closed)/len(core_terms) if core_terms else None
        if label in {'stage1', 'stage2_standard', 'stage2_demotion'} and violations:
            raise RuntimeError(f'{label} left known-taxon violations: {violations[:10]}')
        if before is not None:
            result['before'] = stats(before, tau, complexes)
            result['modified_test_proteins'] = [p for p in annotations if scores.get(p, {}) != before.get(p, {})]
            a, b = active(before, tau), active(scores, tau)
            result['promotions'] = sum(len(b[p] - a.get(p, set())) for p in b)
            result['demotions'] = sum(len(a[p] - b.get(p, set())) for p in a)
        out['runs'].append(result)
    record('before', raw, native)
    native_ann = active(raw, native)
    native_costs, _ = tc.compute_demotion_flip_cost(raw, native)
    out['native_active_terms'] = sorted(set().union(*native_ann.values()))
    out['native_flip_costs'] = dict(native_costs)
    if args.infer_taxons:
        solution, inferred, _ = tc.solve_sat_ortools_taxon_assignment(
            native_ann, in_tax, never_tax, hierarchy, native_costs, tax_isa, tax_disjoint, tax_union)
        unknown = tc.get_final_predictions(solution, raw, hierarchy)
        record('stage1_unknown', unknown, native, raw)
        out['runs'][-1]['solver_taxons'] = inferred
        out['runs'][-1]['retained_active_terms'] = sorted(set().union(*active(unknown, native).values()))
        del unknown
    if args.funtaxis:
        published = json.loads(Path(args.funtaxis).read_text())
        forbidden = go.normalize_terms(published['taxa'][args.taxon]['terms'])
        filtered = {p: {g: s for g, s in scores.items() if g not in forbidden} for p, scores in raw.items()}
        record('published_funtaxis', filtered, native, raw)
        out['funtaxis_manifest_sha256'] = sha(args.funtaxis)
        del filtered
    for tau in ([.001, .01, .05, .1, .2, .3, .4, .5] if args.sweep else [native]):
        if tau != native:
            record('before', raw, tau)
        ann = active(raw, tau)
        costs, _ = tc.compute_demotion_flip_cost(raw, tau)
        start = time.monotonic()
        adjusted, _ = tc.solve_sat_ortools_hierarchy(ann, in_tax, never_tax, hierarchy, costs,
                                                    tax_isa, tax_disjoint, tax_union, args.taxon)
        repaired = tc.get_final_predictions(adjusted, raw, hierarchy)
        record('stage1', repaired, tau, raw)
        out['runs'][-1]['seconds_stage1'] = time.monotonic() - start
        out['runs'][-1]['before'] = stats(raw, tau, complexes)
        if tau == native:
            target = Path(args.output).parent / 'stage1'
            target.mkdir(exist_ok=True)
            cc.save_predictions(repaired, target / pred_file.name.replace('predictions_', 'optimized_'))
        if tau == native and args.aspect == 'cc' and args.stage2:
            for demotion in [False, True]:
                ann2, terms2 = cc.get_annotated_predictions(repaired, tau)
                candidates = cc.get_participating_proteins(ann2, complexes, repaired, top_k=50, demotion_only=demotion)
                pairs, candidates = cc.prepare_sparse_domain(ann2, hierarchy, candidates, complexes)
                costs2, _ = cc.compute_flip_cost(repaired, tau)
                solver_stats = {}
                solution, _ = cc.solve_sat_ortools_optimized(ann2, pairs, hierarchy, costs2, candidates, demotion_only=demotion, time_limit=args.solver_seconds, solver_stats=solver_stats)
                final = cc.get_final_predictions(solution, repaired, terms2, tau)
                label = 'stage2_demotion' if demotion else 'stage2_standard'
                record(label, final, tau, repaired)
                out["runs"][-1]["solver"] = solver_stats
                if out['runs'][-1]['n_singleton_complexes']:
                    raise RuntimeError('Stage 2 left singleton complexes')
                target = Path(args.output).parent / label
                target.mkdir(exist_ok=True)
                cc.save_predictions(final, target / pred_file.name.replace('predictions_', 'optimized_'))
                del final
            if args.stage2_gated:
                final, gated_info = gated_stage2(args, repaired, tau, hierarchy, isa_hierarchy, complexes)
                record('stage2_gated', final, tau, repaired)
                out['runs'][-1]['solver'] = gated_info.pop('solver')
                out['runs'][-1]['gated'] = {k: gated_info[k] for k in ['n_singletons', 'n_with_evidence', 'outcome_counts', 'evidence_inputs']}
                target = Path(args.output).parent / 'stage2_gated'
                target.mkdir(exist_ok=True)
                cc.save_predictions(final, target / pred_file.name.replace('predictions_', 'optimized_'))
                sidecar = Path(args.output).parent / 'gated_variants' / f'{args.taxon}.json.gz'
                sidecar.parent.mkdir(exist_ok=True)
                with gzip.open(sidecar, 'wt') as fh:
                    json.dump(gated_info['sidecar'], fh)
                del final
        del repaired
    if args.soft:
        ann = active(raw, native)
        costs, _ = tc.compute_demotion_flip_cost(raw, native)
        counts = Counter()
        for terms in ann.values():
            counts.update(set().union(*(go.get_ancestors(g) for g in terms)))
        for penalty in [.2, .5, 1.0]:
            violation_costs = {g: penalty * n for g, n in counts.items()}
            solution, _ = tc.solve_sat_ortools_hierarchy(
                ann, in_tax, never_tax, hierarchy, costs, tax_isa, tax_disjoint,
                tax_union, args.taxon, violation_costs=violation_costs)
            soft = tc.get_final_predictions(solution, raw, hierarchy)
            record(f'soft_stage1_{penalty}', soft, native, raw)
            if args.aspect == 'cc' and args.stage2:
                soft_ann, terms = cc.get_annotated_predictions(soft, native)
                candidates = cc.get_participating_proteins(soft_ann, complexes, soft, top_k=50)
                pairs, candidates = cc.prepare_sparse_domain(soft_ann, hierarchy, candidates, complexes)
                costs2, _ = cc.compute_flip_cost(soft, native)
                solver_stats = {}
                solution, _ = cc.solve_sat_ortools_optimized(
                    soft_ann, pairs, hierarchy, costs2, candidates, violation_penalty=penalty, time_limit=args.solver_seconds, solver_stats=solver_stats)
                final = cc.get_final_predictions(solution, soft, terms, native)
                record(f'soft_stage2_{penalty}', final, native, soft)
                out['runs'][-1]['solver'] = solver_stats
                del final
            del soft
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(str(target) + '.part', 'wt') as fh:
        json.dump(out, fh, allow_nan=False)
    Path(str(target) + '.part').replace(target)


GATED_TAU = {'deepgo-se': .3, 'seq_sim_timeset': .001, 'sprof-go': .001}
GATED_PREDICTORS = ['seq_sim_timeset', 'deepgo-se', 'sprof-go']
GENERIC_COMPLEX = {'GO:0032991'}


def _load_selected(path, terms):
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


def gated_stage2(args, scores, tau, H, isa_h, complexes, penalty=100.):
    """Evidence-gated Stage 2 (sx_prefer); mirrors investigate-complex/variants.py."""
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
    others, evidence_inputs = {}, {}
    for op in GATED_PREDICTORS:
        if op == args.predictor:
            continue
        files = sorted((Path(args.gated_evidence.format(predictor=op))).glob(f'optimized_fold_*_taxon_{args.taxon}.tsv'))
        if len(files) != 1:
            raise FileNotFoundError(f'gated evidence {op} {args.taxon}: {files}')
        evidence_inputs[str(files[0])] = sha(files[0])
        sc = _load_selected(files[0], set(singles))
        others[op] = {(g, p) for p, d in sc.items() for g, v in d.items() if v > GATED_TAU[op]}
        del sc
    ev = {}
    for X in singles:
        L = lone[X]
        parents = {y for y in isa_h.get(X, ()) if y in complexes and y not in GENERIC_COMPLEX}
        e = defaultdict(list)
        for p in names:
            if p == L:
                continue
            if parents & ann[p]:
                e[p].append('rel')
            if any((X, p) in o for o in others.values()):
                e[p].append('xpred')
        ev[X] = dict(e)
    cand = {X: sorted((p for p, t in ev[X].items() if 'xpred' in t),
                      key=lambda p: (-len(ev[X][p]), -scores[p].get(X, 0.0), p))[:50] for X in singles}
    demotable = {X for X in singles if not cand[X]}
    costs, _ = cc.compute_flip_cost(scores, tau)
    part = {X: sorted({lone[X]} | set(cand[X])) for X in singles}
    pairs, affected = cc.prepare_sparse_domain(ann, H, part, complexes)
    prom = {(g, q) for X in singles for q in cand[X] for g in {X} | anc(X) if g not in ann[q]}
    dem = {(g, lone[X]) for X in demotable for g in {X} | desc(X) if g in ann[lone[X]]}
    st = {}
    t0 = time.monotonic()
    # As in variants.py (--workers 3): the gated model is solved with 3 CP-SAT
    # workers; with the pipeline's single worker the hard model often stops
    # FEASIBLE at 120 s (e.g. DeepGO-SE 10116: 18.43 vs proven optimum 17.70).
    original_solver = cc.cp_model.CpSolver
    workers = args.gated_workers
    class _Solver(original_solver):
        def Solve(self, model, *x, **k):
            self.parameters.num_search_workers = workers
            return super().Solve(model, *x, **k)
    cc.cp_model.CpSolver = _Solver
    try:
        sol, st = _gated_solve(singles, ann, pairs, H, costs, affected, prom, dem, args, penalty)
    finally:
        cc.cp_model.CpSolver = original_solver
    st['workers'] = workers
    st['wall'] = time.monotonic() - t0
    return _gated_finish(args, scores, terms, tau, ann, complexes, singles, lone, ev, cand, st, sol, evidence_inputs)


def _gated_solve(singles, ann, pairs, H, costs, affected, prom, dem, args, penalty):
    st = {}
    if not singles:
        sol = {}
        st = dict(status='OPTIMAL', seconds=0.0, objective=0.0, bound=0.0, optimal=True, formulation='none')
    else:
        try:
            sol, _ = cc.solve_sat_ortools_optimized(ann, pairs, H, costs, affected, violation_penalty=None,
                                                    time_limit=args.solver_seconds, solver_stats=st,
                                                    allowed_promotions=prom, allowed_demotions=dem)
            st['formulation'] = 'hard'
        except RuntimeError:
            st = {}
            sol, _ = cc.solve_sat_ortools_optimized(ann, pairs, H, costs, affected, violation_penalty=penalty,
                                                    time_limit=args.solver_seconds, solver_stats=st,
                                                    allowed_promotions=prom, allowed_demotions=dem)
            st['formulation'] = 'soft_fallback'
    return sol, st


def _gated_finish(args, scores, terms, tau, ann, complexes, singles, lone, ev, cand, st, sol, evidence_inputs):
    final = cc.get_final_predictions(sol, scores, terms, tau)
    after = active(final, tau)
    promotions = sorted([g, p] for p in after for g in after[p] - ann.get(p, set()))
    demotions = sorted([g, p] for p in ann for g in ann[p] - after.get(p, set()))
    cnt = Counter(g for ts in after.values() for g in ts & complexes)
    outcome = {}
    for X in singles:
        if X not in after.get(lone[X], set()):
            outcome[X] = 'demoted'
        elif cnt[X] >= 2:
            outcome[X] = 'promoted'
        else:
            outcome[X] = 'unresolved'
    variant = dict(solver=dict(st), promotions=promotions, demotions=demotions, outcome=outcome,
                   n_present=sum(1 for g in cnt if cnt[g]), n_singletons_after=sum(1 for g in cnt if cnt[g] == 1))
    sidecar = dict(predictor=args.predictor, taxon=args.taxon, tau=tau, n_singletons=len(singles), dense_floor=None,
                   evidence=ev, lone=lone, lone_score={X: scores[lone[X]].get(X, 0.0) for X in singles},
                   candidates=cand, variants={'sx_prefer': variant})
    return final, dict(solver=st, n_singletons=len(singles), n_with_evidence=sum(bool(cand[X]) for X in singles),
                       outcome_counts=dict(Counter(outcome.values())), evidence_inputs=evidence_inputs, sidecar=sidecar)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--annotations', help='Rebuilt source-reference directory; overrides the legacy input subdirectory')
    for arg in ['data', 'resources', 'predictor', 'aspect', 'taxon', 'output']:
        ap.add_argument('--' + arg, required=True)
    ap.add_argument('--proteins')
    ap.add_argument('--funtaxis')
    ap.add_argument('--sweep', action='store_true')
    ap.add_argument('--stage2', action='store_true')
    ap.add_argument('--soft', action='store_true')
    ap.add_argument('--infer-taxons', action='store_true')
    ap.add_argument('--solver-seconds', type=float, default=120)
    ap.add_argument('--gated-workers', type=int, default=3)
    ap.add_argument('--stage2-gated', action='store_true', help='also run the evidence-gated Stage 2 (sx_prefer)')
    ap.add_argument('--gated-evidence', default=None,
                    help='directory template with {predictor} holding other predictors\' Stage 1 optimized_*.tsv')
    run(ap.parse_args())
