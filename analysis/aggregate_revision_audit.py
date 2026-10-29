"""Aggregate a complete task manifest; never silently omit a failed organism."""
import argparse
from collections import defaultdict, Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaef_patches.ontology import Ontology, NAMESPACES, FUNC_DICT
from pipeline.genome_scores_evaluator import evaluate_prediction_metrics, load_annotations
from analysis.evaluation_scope import load_temporal_proteins, select_annotations


def score_delta(base, scores):
    changes = {g: value for g, value in scores.items() if g not in base or value != base[g]}
    changes.update({g: None for g in base if g not in scores})
    return changes


def apply_score_delta(base, delta):
    result = {g: value for g, value in base.items() if delta.get(g, value) is not None}
    result.update({g: value for g, value in delta.items() if value is not None})
    return result


def checked_reference_labels(data, directory, reference, temporal):
    """Reuse label-independent repairs only with complete stored score coverage."""
    taxon = data['taxon']
    filename = reference['reference_file']
    captured = [sha for name, sha in data['inputs'].items() if Path(name).name == filename]
    if len(captured) != 1 or captured[0] not in {reference['original_sha256'], reference['reference_sha256']}:
        raise ValueError(f'Unidentified original/reference annotation input for {taxon}')
    mask = None if data['predictor'].endswith('_heldout') else temporal
    labels = select_annotations(load_annotations(Path(directory)/filename), taxon, mask)
    for run in data['runs']:
        missing = set(labels)-set(run['test_scores'])
        if missing:
            raise ValueError(f'Rerun required: {data["predictor"]}/{data["aspect"]}/{taxon} lacks scores for {sorted(missing)[:10]}')
    return labels


def aggregate(args):
    tasks = [tuple(line.split()) for line in Path(args.tasks).read_text().splitlines() if line.strip()]
    if args.predictor:
        tasks = [t for t in tasks if t[0] == args.predictor]
    if args.aspect:
        tasks = [t for t in tasks if t[1] == args.aspect]
    if not tasks or any(len(t) != 3 for t in tasks):
        raise ValueError('Empty or malformed task manifest')
    if len(tasks) != len(set(tasks)):
        raise ValueError('Duplicate task in manifest')
    missing = [str(Path(args.results) / model / aspect / f'{taxon}.json.gz')
               for model, aspect, taxon in tasks
               if not (Path(args.results) / model / aspect / f'{taxon}.json.gz').exists()]
    if missing:
        raise FileNotFoundError(f'{len(missing)} required task outputs missing; first: {missing[:5]}')
    go = Ontology(args.go)
    temporal = load_temporal_proteins(args.proteins)
    heldout_train = set(Path(args.heldout_train).read_text().split())
    counts = {'timeset': Counter(), 'heldout': Counter()}
    ancestor_cache = {}
    ic_reference_inputs = {}
    reference_path = Path(args.annotations)/'reference-manifest.json'
    if not reference_path.is_file():
        raise FileNotFoundError('Reference labels must be rebuilt from DAT sources; reference-manifest.json is required')
    reference_manifest = json.loads(reference_path.read_text())
    if reference_manifest['go_sha256'] != hashlib.sha256(Path(args.go).read_bytes()).hexdigest():
        raise ValueError('Reference label build used a different GO ontology')
    reference_records = {r['taxon']:r for r in reference_manifest['records']}
    expected_files = {r['reference_file']:r['reference_sha256'] for r in reference_manifest['records']}
    if len(reference_records) != len(reference_manifest['records']) or len(expected_files) != len(reference_records):
        raise ValueError('Duplicate reference input records')
    def ancestors(term):
        if term not in ancestor_cache:
            ancestor_cache[term] = go.get_ancestors(term)
        return ancestor_cache[term]
    for path in sorted(Path(args.annotations).glob('annots_taxon_*.tsv')):
        ic_reference_inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected_files.get(path.name) != ic_reference_inputs[str(path)]:
            raise ValueError(f'Reference input checksum differs: {path}')
        taxon = path.stem.split('_')[-1]
        for protein, terms in load_annotations(path).items():
            closed = set().union(*(ancestors(g) for g in terms))
            if protein not in temporal.get(taxon, ()):
                counts['timeset'].update(closed)
            if taxon in heldout_train:
                counts['heldout'].update(closed)
    if {Path(p).name for p in ic_reference_inputs} != set(expected_files):
        raise ValueError('Reference input manifest is incomplete')
    ic = {}
    for split, corpus in counts.items():
        go.calculate_ic_from_counts(corpus)
        ic[split] = (go.ic.copy(), go.ic_norm)
    groups = defaultdict(list)
    baseline_scores = {}
    annotations = defaultdict(dict)
    source_records = []
    code_manifest = None
    scoped_counts = defaultdict(Counter)
    identifier_counts = defaultdict(Counter)
    for model, aspect, taxon in tasks:
        path = Path(args.results) / model / aspect / f'{taxon}.json.gz'
        with gzip.open(path, 'rt') as fh:
            data = json.load(fh)
        if (data['predictor'], data['aspect'], data['taxon']) != (model, aspect, taxon):
            raise ValueError(f'Task identity mismatch: {path}')
        if data.get('schema_version') != 2:
            raise ValueError(f'Unsupported result schema: {path}')
        if code_manifest is None:
            code_manifest = data['code']
        elif data['code'] != code_manifest:
            raise ValueError(f'Mixed code versions: {path}')
        config = data['configuration']
        if config.get('input_alignment') != 'primary_ids_max_score_is_a_part_of_closure_v1':
            raise ValueError(f'Incompatible input-alignment convention: {path}')
        if 'identifier_normalization' not in data:
            raise ValueError(f'Result predates canonical GO identifier handling: {path}')
        identifier_counts[model, aspect].update(data['identifier_normalization'])
        native = data['native_tau']
        taus = [.001, .01, .05, .1, .2, .3, .4, .5] if config['sweep'] else [native]
        expected_modes = {(mode, tau) for mode in ['before', 'stage1'] for tau in taus}
        if config['infer_taxons']:
            expected_modes.add(('stage1_unknown', native))
        if 'funtaxis_manifest_sha256' in data:
            expected_modes.add(('published_funtaxis', native))
        if config['stage2'] and aspect == 'cc':
            expected_modes.update((mode, native) for mode in ['stage2_standard', 'stage2_demotion'])
        if config['soft']:
            for penalty in [.2, .5, 1.0]:
                expected_modes.add((f'soft_stage1_{penalty}', native))
                if config['stage2'] and aspect == 'cc':
                    expected_modes.add((f'soft_stage2_{penalty}', native))
        if {(r['mode'], r['tau']) for r in data['runs']} != expected_modes:
            raise ValueError(f'Missing or unexpected conditions: {path}')
        labels = checked_reference_labels(data, args.annotations, reference_records[taxon], temporal)
        scoped_counts[model, aspect].update({
            'n_mask_proteins':len(labels) if model.endswith('_heldout') else len(temporal.get(taxon, ())),
            'n_annotation_matches':len(labels)})
        annotations[model, aspect].update({taxon + ':' + p: set(ts) for p, ts in labels.items()})
        baseline = next(r['test_scores'] for r in data['runs'] if r['mode'] == 'before' and r['tau'] == native)
        baseline_scores[model, aspect, taxon] = baseline
        seen = set()
        for run in data['runs']:
            scores = run.pop('test_scores')
            run['test_score_deltas'] = {p: delta for p, terms in scores.items() if (delta := score_delta(baseline.get(p, {}), terms))}
            key = (model, aspect, run['mode'], run['tau'])
            if key in seen:
                raise ValueError(f'Duplicate run in {path}: {key}')
            seen.add(key)
            groups[key].append((taxon, run))
        source_records.append({'task': [model, aspect, taxon], 'inputs': data['inputs'],
                               'result_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                               'code': data.get('code'), 'software': data['software'], 'configuration': data['configuration'],
                               'protein_scope': data['protein_scope']})
    rows = []
    baseline_metrics = {}
    for (model, aspect, mode, tau), records in sorted(groups.items()):
        expected = sum(m == model and a == aspect for m, a, _ in tasks)
        if len(records) != expected:
            raise ValueError(f'Incomplete condition {model}/{aspect}/{mode}/{tau}: {len(records)}/{expected}')
        split = 'heldout' if model.endswith('_heldout') else 'timeset'
        go.ic, go.ic_norm = ic[split]
        predictions = {}
        for taxon, run in records:
            for protein, base in baseline_scores[model, aspect, taxon].items():
                scores = apply_score_delta(base, run['test_score_deltas'].get(protein, {}))
                propagated = defaultdict(float)
                for term, score in scores.items():
                    for ancestor in ancestors(term):
                        propagated[ancestor] = max(propagated[ancestor], score)
                predictions[taxon + ':' + protein] = dict(propagated)
        ann = annotations[model, aspect]
        # The repair threshold changes structural counts, not raw continuous
        # scores. Verify that invariant before reusing their expensive PR curve.
        if mode == 'before':
            if any(run['test_score_deltas'] for _, run in records):
                raise ValueError('Raw baseline scores differ across repair thresholds')
            if (model, aspect) not in baseline_metrics:
                baseline_metrics[model, aspect] = evaluate_prediction_metrics(
                    go, aspect, ann, predictions, reuse_ic=True)
            metrics = baseline_metrics[model, aspect]
        else:
            metrics = evaluate_prediction_metrics(go, aspect, ann, predictions, reuse_ic=True)
        namespace = go.get_namespace_terms(NAMESPACES[aspect]) - {FUNC_DICT[aspect]}
        valid = {p for p, ts in ann.items() if set().union(*(ancestors(g) for g in ts)) & namespace}
        changed = {taxon+':'+p for taxon, r in records for p in r.get('modified_test_proteins', ())} & valid
        modified_metrics = evaluate_prediction_metrics(go, aspect, {p: ann[p] for p in changed}, predictions, reuse_ic=True) if changed else None
        if metrics is None:
            raise ValueError(f'No assessable metrics for {model}/{aspect}/{mode}')
        active = sum(r.get('before', r)['active_predictions'] for _, r in records)
        grid = sum(r.get('before', r)['grid_decisions'] for _, r in records)
        demotions = sum(r.get('demotions', 0) for _, r in records)
        promotions = sum(r.get('promotions', 0) for _, r in records)
        row = {'predictor': model, 'aspect': aspect, 'mode': mode, 'tau': tau,
               'n_organisms': expected, 'n_annotated_proteins': len(valid),
               'n_annotated_organisms': len({p.split(':',1)[0] for p in valid}),
               'scope_counts': dict(scoped_counts[model, aspect]),
               'identifier_normalization': dict(identifier_counts[model, aspect]),
               'n_modified_test_proteins': len(changed), 'modified_metrics': modified_metrics,
               'active_input': active, 'grid_input': grid, 'demotions': demotions, 'promotions': promotions,
               'demotion_active_pct': 100 * demotions / active if active else None,
               'flip_grid_pct': 100 * (demotions + promotions) / grid if grid else None,
               'known_taxon_consistency': sum(r['known_taxon_consistent'] for _, r in records) / expected,
               'complex_coherence': sum(r['complex_coherence'] for _, r in records) / expected if aspect == 'cc' else None,
               'gaef_taxon_consistency': sum(r['gaef_taxon_consistent'] for _, r in records) / expected,
               'solver_status_counts': dict(Counter(r['solver']['status'] for _, r in records if 'solver' in r)),
               'solver_objective_sum': sum(r['solver']['objective'] for _, r in records if 'solver' in r),
               'solver_bound_sum': sum(r['solver']['bound'] for _, r in records if 'solver' in r),
               'metrics': metrics}
        for key in ['process_coherence', 'pathway_coherence', 'essential_completeness']:
            values = [r[key] for _, r in records if r.get(key) is not None]
            row[key] = sum(values)/len(values) if values else None
            row[key+'_n'] = len(values)
        rows.append(row)
    Path(args.output).write_text(json.dumps({'schema_version': 2,
        'ic_definition': 'log2(min parent count / term count), fixed reference: timeset excludes entry-date test proteins; heldout uses training taxa only',
        'sources': source_records, 'rows': rows,
        'ic_reference_inputs': ic_reference_inputs,
        'reference_manifest_sha256': hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        'reference_label_policy': 'Source-DAT evidence with is_a/part_of closure; every evaluation protein has a stored score vector in every condition',
        'ic_split_inputs': {str(Path(p)): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                            for p in [args.proteins, args.heldout_train, args.go]},
        'aggregation_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['tasks', 'results', 'go', 'annotations', 'proteins', 'output']:
        ap.add_argument('--' + name, required=True)
    ap.add_argument('--heldout-train', required=True)
    ap.add_argument('--predictor')
    ap.add_argument('--aspect', choices=['cc', 'mf', 'bp'])
    aggregate(ap.parse_args())
