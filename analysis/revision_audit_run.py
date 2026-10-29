"""Rerun one predictor/aspect/organism with explicit inputs and test membership.

Each output is independent and carries input checksums. Aggregation is a
separate step, so missing organisms cannot silently produce a complete table.
"""
import argparse
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

ROOT = Path(__file__).resolve().parents[1]
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
                             'input_alignment': 'primary_ids_max_score_is_a_part_of_closure_v1'},
           'software': {'python': platform.python_version(), **{pkg: importlib.metadata.version(pkg) for pkg in ['ortools', 'numpy', 'pandas', 'scikit-learn']}}, 'predictor': args.predictor, 'aspect': args.aspect,
           'taxon': args.taxon, 'native_tau': native, 'n_complex_terms': len(complexes),
           'protein_scope': scope_manifest(args.proteins) if args.proteins else {'definition': 'held-out organism'},
           'inputs': {str(f): sha(f) for f in [pred_file, ann_file, go_file, constraints_file, ROOT / 'complex_coherence/protein_complexes.tsv'] + tax_files + sorted(cdir.glob('*')) if f.is_file()},
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
    run(ap.parse_args())
