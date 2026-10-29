"""Evaluate seven-rank lineage inference from complete per-task rerun outputs.

Taxa forced by retained GO terms are distinguished from arbitrary solver
assignments. Forcedness is conditional on the selected minimum-cost GO repair;
alternative equally optimal GO repairs are not enumerated.
"""
import argparse
import csv
import hashlib
from collections import defaultdict, Counter
import gzip
import json
from pathlib import Path
import sys
from ortools.sat.python import cp_model

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from taxon_consistency import adjust_ortools as tc
from gaef_patches.ontology import Ontology
from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage, evaluate_single_lineage, RANKS_7


def supported(terms, only, never, isa, disjoint, unions, go, ranks, forced=True):
    closed = set().union(*(go.get_ancestors(g) for g in terms))
    yes = {t for g in closed for t in only.get(g, ())}
    no = {t for g in closed for t in never.get(g, ())}
    model = cp_model.CpModel()
    variables = tc._add_taxon_model(model, yes | no, isa, disjoint, unions)
    for taxon in yes:
        model.Add(variables[taxon] == 1)
    for taxon in no:
        model.Add(variables[taxon] == 0)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 42
    status = solver.Solve(model)
    if status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        if status != cp_model.INFEASIBLE:
            raise RuntimeError(f'Undecided lineage: {solver.StatusName(status)}')
        return None
    if not forced:
        return True
    required = set()
    # Only ranked nodes affect the requested evaluation. Check logical
    # entailment, rather than treating a chosen union member as uniquely forced.
    for taxon, variable in variables.items():
        if ranks.get(taxon) not in RANKS_7:
            continue
        model.ClearAssumptions()
        model.AddAssumption(variable.Not())
        status = solver.Solve(model)
        if status == cp_model.INFEASIBLE:
            required.add(taxon)
        elif status not in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
            raise RuntimeError('Undecided forced taxon')
    return required


def run(args):
    resources = Path(args.resources)
    parents, ranks, names = parse_ncbitaxon_obo(resources/'ncbitaxon.obo')
    only, never = tc.load_constraints(resources/'go_taxon_constraints_extracted_obo.tsv')
    isa, disjoint, unions = tc.load_taxon_hierarchy(resources/'taxon_hierarchy.tsv', resources/'ncbitaxon_hierarchy.tsv')
    go = Ontology(resources/'go-basic.obo')
    only, never = go.normalize_constraints(only), go.normalize_constraints(never)
    hierarchy = {g: set(obj['is_a']) for g, obj in go.ont.items()}
    timeset_taxa = sorted(Path(args.taxa).read_text().split())
    heldout_taxa = sorted(Path(args.heldout_taxa).read_text().split())
    predictors = ['deepgo-se', 'seq_sim_timeset', 'sprof-go', 'seq_sim_heldout', 'mlp_heldout']
    records = []
    input_hashes = {}
    for predictor in predictors:
        taxa = heldout_taxa if predictor.endswith('_heldout') else timeset_taxa
        for taxon in taxa:
            true = get_true_lineage(taxon, parents)
            terms_by_aspect, native, costs = {}, {}, defaultdict(float)
            aspect_evals = {}
            for aspect in ['cc', 'mf', 'bp']:
                path = Path(args.results)/predictor/aspect/f'{taxon}.json.gz'
                with gzip.open(path, 'rt') as fh:
                    data = json.load(fh)
                input_hashes[f'{predictor}/{aspect}/{taxon}'] = hashlib.sha256(path.read_bytes()).hexdigest()
                result = next(r for r in data['runs'] if r['mode'] == 'stage1_unknown')
                terms_by_aspect[aspect] = set(result['retained_active_terms'])
                native[aspect] = set(data['native_active_terms'])
                for g, cost in data['native_flip_costs'].items():
                    costs[g] += cost
                required = supported(terms_by_aspect[aspect], only, never, isa, disjoint, unions, go, ranks)
                if required is None:
                    raise ValueError(f'Inconsistent unknown repair: {path}')
                evaluation = evaluate_single_lineage(true, required, ranks, parents)
                aspect_evals[aspect] = evaluation
                records.append({'predictor': predictor, 'taxon': taxon, 'aspect': aspect,
                    'supported': evaluation,
                    'solver_assignment': evaluate_single_lineage(true, set(result['solver_taxons']), ranks, parents)})
            joint, inferred, _ = tc.solve_sat_ortools_taxon_assignment(native, only, never, hierarchy, costs, isa, disjoint, unions)
            retained = {g for assignments in joint.values() for g, keep in assignments.items() if keep}
            required = supported(retained, only, never, isa, disjoint, unions, go, ranks)
            if required is None:
                raise ValueError('Inconsistent joint repair')
            rank_agreement = {}
            for rank in RANKS_7:
                observed = [set(ev['inferred_ranks'].get(rank, ())) for ev in aspect_evals.values()]
                nonempty = [x for x in observed if x]
                rank_agreement[rank] = {'resolved_aspects': len(nonempty),
                    'compatible': len(set().union(*nonempty)) <= 1,
                    'all_three_resolved_equal': len(nonempty) == 3 and all(x == nonempty[0] for x in nonempty)}
            records.append({'predictor': predictor, 'taxon': taxon, 'aspect': 'joint',
                'supported': evaluate_single_lineage(true, required, ranks, parents),
                'solver_assignment': evaluate_single_lineage(true, set(inferred), ranks, parents),
                'cross_aspect_rank_agreement': rank_agreement,
                'cross_aspect_constraints_compatible': bool(supported(set().union(*terms_by_aspect.values()), only, never, isa, disjoint, unions, go, ranks, forced=False))})
    summaries = []
    for predictor in predictors:
        for aspect in ['cc', 'mf', 'bp', 'joint']:
            subset = [r for r in records if r['predictor'] == predictor and r['aspect'] == aspect]
            row = {'predictor': predictor, 'aspect': aspect, 'n': len(subset)}
            for source in ['supported', 'solver_assignment']:
                row[source] = {key: dict(Counter(r[source][key] for r in subset)) for key in ['deepest_correct_rank', 'first_divergence_rank']}
                row[source]['exact_prefix_count'] = sum(r[source]['is_exact_prefix'] for r in subset)
                row[source]['unresolved_count'] = sum(r[source]['n_inferred_nodes'] == 0 for r in subset)
            if aspect == 'joint':
                row['cross_aspect_constraints_compatible'] = sum(r['cross_aspect_constraints_compatible'] for r in subset)
                row['cross_aspect_rank_agreement'] = {rank: {key: sum(r['cross_aspect_rank_agreement'][rank][key] for r in subset) for key in ['compatible', 'all_three_resolved_equal']} for rank in RANKS_7}
            summaries.append(row)
    with Path(args.output).with_suffix('.tsv').open('w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['predictor', 'taxon', 'aspect', 'rank', 'true_id', 'true_name', 'required_ids', 'required_names', 'solver_ids'])
        for record in records:
            for rank in RANKS_7:
                truth = record['supported']['true_ranks'].get(rank, '')
                inferred = record['supported']['inferred_ranks'].get(rank, ())
                writer.writerow([record['predictor'], record['taxon'], record['aspect'], rank, truth, names.get(truth, ''),
                    ';'.join(inferred), ';'.join(names.get(t, t) for t in inferred), ';'.join(record['solver_assignment']['inferred_ranks'].get(rank, ()))])
    resource_hashes = {name: hashlib.sha256((resources/name).read_bytes()).hexdigest()
        for name in ['ncbitaxon.obo', 'go-basic.obo', 'taxon_hierarchy.tsv',
                     'ncbitaxon_hierarchy.tsv', 'go_taxon_constraints_extracted_obo.tsv']}
    Path(args.output).write_text(json.dumps({'analysis_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256': input_hashes, 'resource_sha256': resource_hashes,
        'tie_policy': __doc__, 'rank_order': RANKS_7, 'summary': summaries, 'records': records}, indent=2)+'\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['results', 'resources', 'taxa', 'heldout-taxa', 'output']:
        ap.add_argument('--'+name, required=True)
    run(ap.parse_args())
