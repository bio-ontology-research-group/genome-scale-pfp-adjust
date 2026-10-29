"""Reproducible extension experiments; see reproducibility/EXTENSIONS.md."""
import argparse
from collections import Counter, defaultdict
import contextlib
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from complex_coherence import adjust_ortools as cc
from complex_coherence.conservative import repair
from taxon_consistency import adjust_ortools as tc
from gaef_patches.ontology import Ontology
from analysis.audit_lineage_results import supported
from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage, evaluate_single_lineage, RANKS_7
from analysis.complex_validation.build_subunit_index import parse_complexportal_file
from analysis.complex_validation.validate_stage2 import protein_to_accession


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def resources(path):
    go = Ontology(path/'go-basic.obo')
    hierarchy = {g: set(o['is_a']) for g, o in go.ont.items()}
    exclusions = go.normalize_terms(cc.load_homodimer_terms(ROOT/'complex_coherence/protein_complexes.tsv'))
    isa_go = Ontology(path/'go-basic.obo', with_rels=False)
    isa_hierarchy = {g: set(o['is_a']) for g, o in isa_go.ont.items()}
    complexes = cc.get_heteromeric_complexes(exclusions, set(), isa_hierarchy)
    return go, hierarchy, complexes


def lineage(args):
    go, hierarchy, _ = resources(args.resources)
    only, never = tc.load_constraints(args.resources/'go_taxon_constraints_extracted_obo.tsv')
    only, never = go.normalize_constraints(only), go.normalize_constraints(never)
    isa, disjoint, unions = tc.load_taxon_hierarchy(args.resources/'taxon_hierarchy.tsv', args.resources/'ncbitaxon_hierarchy.tsv')
    parents, ranks, _ = parse_ncbitaxon_obo(args.resources/'ncbitaxon.obo')
    query = {t for t, rank in ranks.items() if rank in RANKS_7}
    native, costs, inputs, previous = {}, {}, {}, {}
    for aspect in ['cc','mf','bp']:
        path = args.results/args.predictor/aspect/f'{args.taxon}.json.gz'
        inputs[str(path)] = sha(path)
        with gzip.open(path, 'rt') as f:
            data = json.load(f)
        native[aspect] = set(data['native_active_terms'])
        costs[aspect] = data['native_flip_costs']
        previous[aspect] = set(next(r for r in data['runs'] if r['mode']=='stage1_unknown')['retained_active_terms'])
        del data
    joint_costs = defaultdict(float)
    for c in costs.values():
        for g, v in c.items():
            joint_costs[g] += v
    records = []
    truth = get_true_lineage(args.taxon, parents)
    for aspect in ['cc','mf','bp','joint']:
        annotations = native if aspect == 'joint' else {aspect: native[aspect]}
        weight = joint_costs if aspect == 'joint' else costs[aspect]
        audit = {}
        result, assignment, _ = tc.solve_sat_ortools_taxon_assignment(
            annotations, only, never, hierarchy, weight, isa, disjoint, unions,
            uncertainty=audit, query_taxa=query, time_limit=args.time_limit)
        retained = {g for terms in result.values() for g, keep in terms.items() if keep}
        conditional = supported(retained, only, never, isa, disjoint, unions, go, ranks)
        required = set(audit['required'])
        if not required <= conditional:
            raise AssertionError('Robust taxa must be required by the selected repair')
        previous_required = (supported(previous[aspect], only, never, isa, disjoint, unions, go, ranks)
                             if aspect != 'joint' else conditional)
        if not required <= previous_required:
            raise AssertionError('Robust taxa disagree with the previously validated optimal repair')
        # A smaller result artifact; absent taxonomy nodes are not assertions.
        outside_count = len(audit.pop('outside_model'))
        records.append(dict(aspect=aspect, uncertainty=audit,
                            taxonomy_nodes_outside_model=outside_count,
                            conditional_required=sorted(conditional),
                            previous_required=sorted(previous_required),
                            removed_conditional=sorted(conditional-required),
                            conditional_evaluation=evaluate_single_lineage(truth, conditional, ranks, parents),
                            robust_evaluation=evaluate_single_lineage(truth, required, ranks, parents)))
    return dict(predictor=args.predictor, taxon=args.taxon, inputs=inputs, records=records)


def standard(scores, hierarchy, complexes, tau, mode, time_limit):
    active, terms = cc.get_annotated_predictions(scores, tau)
    part = cc.get_participating_proteins(active, complexes, scores, top_k=50,
                                        demotion_only=mode=='demotion_only')
    pairs, affected = cc.prepare_sparse_domain(active, hierarchy, part, complexes)
    costs, _ = cc.compute_flip_cost(scores, tau)
    stats = {}
    solution, _ = cc.solve_sat_ortools_optimized(active, pairs, hierarchy, costs, affected,
                        demotion_only=mode=='demotion_only', time_limit=time_limit, solver_stats=stats)
    return cc.get_final_predictions(solution, scores, terms, tau), stats


def load_selected(path, terms):
    """Stream a dense score file, retaining the complete protein universe."""
    scores = {}
    with open(path) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            p = cols[0]
            if p in scores:
                raise ValueError(f'Duplicate protein {p}')
            values = {}
            for token in cols[1:]:
                if '|' in token:
                    g, s = token.split('|')
                    if g in terms:
                        values[g] = float(s)
            scores[p] = values
    return scores


def recovery_case(original, term, seed, tau, hierarchy, complexes, time_limit, corroborated=None):
    """The hidden member is intentionally not an argument to repair inputs.

    Retain one supplied positive. All other target scores are monotonically
    compressed into [0,tau], creating a singleton challenge. Original scores
    provide candidate evidence, never the masked curated membership.
    """
    ancestors = {term} | cc.get_ancestors(term, hierarchy)
    local_hierarchy = {g: hierarchy.get(g, set()) & ancestors for g in ancestors}
    local_complexes = complexes & ancestors
    scores = {p: {g: (max(.9, tau+.001) if p == seed else tau*s.get(term, 0))
                  for g in ancestors} for p, s in original.items()}
    evidence = {(term,p): 'original predictor score above native threshold'
                for p, s in original.items() if p != seed and s.get(term, 0) > tau}
    outputs = {}
    modes = ['standard','demotion_only','conservative']
    if corroborated is not None:
        modes.append('conservative_corroborated')
    for mode in modes:
        if mode.startswith('conservative'):
            eligible = evidence if mode == 'conservative' else {pair:source for pair,source in corroborated.items() if pair[0]==term and pair[1]!=seed}
            final, audit = repair(scores, local_hierarchy, local_complexes, eligible,
                protected_pairs={(term,seed)}, threshold=tau, time_limit=time_limit)
            stats = audit['solver']
        else:
            final, stats = standard(scores, local_hierarchy, local_complexes, tau, mode, time_limit)
        members = sorted(p for p, s in final.items() if s.get(term, 0)>tau)
        outputs[mode] = dict(members=members, seed_removed=seed not in members,
                             unresolved=len(members)==1, solver=stats)
    ranked = sorted((p for g,p in evidence), key=lambda p: (-original[p].get(term, 0),p))
    outputs['highest_score'] = dict(members=[seed]+ranked[:1], seed_removed=False,
                                    unresolved=not ranked)
    if corroborated is not None:
        ranked_cross = sorted((p for g,p in corroborated if g==term and p!=seed),
                              key=lambda p: (-original[p].get(term, 0),p))
        outputs['highest_score_corroborated'] = dict(members=[seed]+ranked_cross[:1],
                                                    seed_removed=False, unresolved=not ranked_cross)
    return outputs, sorted(p for g,p in evidence)


def recovery(args):
    go, hierarchy, complexes = resources(args.resources)
    paths = sorted((args.results/args.predictor/'cc'/'stage1').glob(f'*_taxon_{args.taxon}.tsv'))
    if len(paths)!=1:
        raise ValueError('Expected one Stage 1 file')
    path = paths[0]
    mapping = json.loads(args.id_map.read_text())
    cp_path = args.complexportal/f'{args.taxon}.tsv'
    if not cp_path.exists():
        return dict(predictor=args.predictor, taxon=args.taxon, records=[], exclusion='no_complexportal_taxon')
    candidates = set()
    excluded = Counter()
    for taxon,g,subunits in parse_complexportal_file(cp_path):
        term = go.canonical_id(g)
        if taxon != args.taxon or term not in complexes:
            excluded['outside_operational_vocabulary'] += 1
        elif len(subunits)!=2:
            excluded['not_two_distinct_accessions'] += 1
        else:
            candidates.add((term, tuple(sorted(subunits))))
    original = load_selected(path, {g for g, _ in candidates})
    inverse = defaultdict(list)
    for p in original:
        accession = protein_to_accession(p, mapping)
        if accession:
            inverse[accession].append(p)
    tau = .3 if args.predictor=='deepgo-se' else .001
    cross_evidence, supporting_inputs = corroborating_evidence(args, {g for g,_ in candidates})
    cross_evidence = {pair:source for pair,source in cross_evidence.items() if pair[1] in original}
    records = []
    for term, accessions in sorted(candidates):
        if any(len(inverse[a])!=1 for a in accessions):
            excluded['missing_or_ambiguous_protein_mapping'] += 1
            continue
        proteins = [inverse[a][0] for a in accessions]
        for seed, hidden in [proteins, proteins[::-1]]:
            outputs, eligible = recovery_case(original, term, seed, tau, hierarchy, complexes, args.time_limit, cross_evidence)
            for output in outputs.values():
                output['hidden_recovered'] = hidden in output['members']
                output['other_proposals'] = sorted(set(output['members'])-{seed,hidden})
            records.append(dict(term=term, seed=seed, hidden=hidden,
                                original_hidden_score=original[hidden].get(term,0),
                                hidden_has_evidence=hidden in eligible,
                                hidden_has_corroboration=(term,hidden) in cross_evidence,
                                eligible_candidates=len(eligible), outputs=outputs))
    return dict(predictor=args.predictor, taxon=args.taxon, threshold=tau,
                n_proteins=len(original), candidate_complexes=len(candidates),
                excluded=dict(excluded), records=records, supporting_inputs=supporting_inputs,
                inputs={str(p): sha(p) for p in [path, cp_path, args.id_map]})


def corroborating_evidence(args, singleton_terms):
    """Require agreement of both other predictors; never consult curated labels."""
    predictors = {'deepgo-se','seq_sim_timeset','sprof-go'} - {args.predictor}
    supported_pairs = []
    sources = {}
    for predictor in sorted(predictors):
        paths = sorted((args.results/predictor/'cc'/'stage1').glob(f'*_taxon_{args.taxon}.tsv'))
        if len(paths) != 1:
            raise ValueError(f'Expected one corroborating prediction file for {predictor}')
        scores = load_selected(paths[0], singleton_terms)
        tau = .3 if predictor == 'deepgo-se' else .001
        supported_pairs.append({(g,p) for p, values in scores.items() for g,v in values.items() if v > tau})
        sources[str(paths[0])] = sha(paths[0])
    return {pair: 'agreement of both other predictors at native thresholds'
            for pair in set.intersection(*supported_pairs)}, sources


def conservative(args):
    go, hierarchy, complexes = resources(args.resources)
    paths = sorted((args.results/args.predictor/'cc'/'stage1').glob(f'*_taxon_{args.taxon}.tsv'))
    if len(paths)!=1:
        raise ValueError('Expected one Stage 1 file')
    path = paths[0]
    scores = cc.load_predictions(path)
    tau = .3 if args.predictor=='deepgo-se' else .001
    active, _ = cc.get_annotated_predictions(scores, tau)
    # This experiment assesses a transparent score-based gate. Nonzero score
    # support is not independent experimental evidence and is labelled as such.
    counts = Counter(g for terms in active.values() for g in terms & complexes)
    singletons = {g for g,n in counts.items() if n == 1}
    evidence = {(g,p): 'nonzero upstream score'
                for p,s in scores.items() for g,v in s.items()
                if g in singletons and 0 < v <= tau}
    supporting_inputs = {}
    if args.mode == 'corroborated':
        evidence, supporting_inputs = corroborating_evidence(args, singletons)
        evidence = {pair:source for pair,source in evidence.items()
                    if pair[1] in scores and pair[0] not in active[pair[1]]}
    final, audit = repair(scores, hierarchy, complexes, evidence,
                         threshold=tau, time_limit=args.time_limit)
    dest = args.output.with_suffix('.tsv')
    cc.save_predictions(final, dest)
    audit.update(predictor=args.predictor, taxon=args.taxon,
                 input_sha256=sha(path), export_sha256=sha(dest),
                 initial_singletons=sum(sum(g in terms for terms in active.values())==1 for g in complexes),
                 evidence_pairs=len(evidence), supporting_inputs=supporting_inputs,
                 evidence_policy=('agreement_of_both_other_predictors; no_demotions' if args.mode == 'corroborated'
                                  else 'nonzero_upstream_score; no_demotions'))
    return audit


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode', choices=['lineage','recovery','conservative','corroborated'])
    for name in ['resources','results','output','id-map','complexportal']:
        ap.add_argument('--'+name, type=Path, required=name in ['resources','results','output'])
    ap.add_argument('--predictor', required=True)
    ap.add_argument('--taxon', required=True)
    ap.add_argument('--time-limit', type=float, default=120)
    args = ap.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with contextlib.redirect_stdout(sys.stderr):
        result = (conservative if args.mode == 'corroborated' else globals()[args.mode])(args)
    result['source_sha256'] = {str(p.resolve().relative_to(ROOT)):sha(p) for p in
            [Path(__file__), ROOT/'complex_coherence/conservative.py',
             ROOT/'complex_coherence/adjust_ortools.py',ROOT/'taxon_consistency/adjust_ortools.py']}
    result['resource_sha256'] = {str(p):sha(p) for p in args.resources.glob('*')
                               if p.name in ['go-basic.obo','ncbitaxon.obo','go_taxon_constraints_extracted_obo.tsv',
                                             'taxon_hierarchy.tsv','ncbitaxon_hierarchy.tsv']}
    temporary = args.output.with_suffix('.partial')
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    temporary.replace(args.output)

if __name__=='__main__':
    main()
