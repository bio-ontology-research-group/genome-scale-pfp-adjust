"""Evaluate a threshold-independent known-taxon filter with the paper's evaluator.

Reads the accepted corrected per-task outputs ($ACCEPTED/<pred>/<aspect>/<taxon>.json.gz),
rebuilds labels, IC corpus and protein scope EXACTLY as analysis/aggregate_revision_audit.py
does, and evaluates (pipeline.genome_scores_evaluator.evaluate_prediction_metrics):

  before          raw canonical, propagated scores (accepted 'before' run, native tau)
  stage1          accepted Stage 1 export (known taxon, native tau)
  stage1_unknown  accepted Stage 1 export with inferred lineage (if present)
  filter_all      'before' with EVERY known-taxon-impossible term removed (score 0) at every
                  score, i.e. any term whose is_a/part_of closure contains a term whose
                  only_in_taxon / never_in_taxon constraint excludes the organism's lineage.
                  The violation test is the same one revision_audit_run.py uses to certify
                  Stage 1 ('known_taxon_consistent').

Metrics are reported on all assessable proteins and on the proteins whose scores the
condition changed. For each evaluation we also capture the 101 macro (precision, recall)
points and recompute AUPR with the ORIGINAL submission's rule (np.trapz over
argsort(recall), no (0,1) endpoint) for attribution.
"""
import argparse
from collections import defaultdict, Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ap = argparse.ArgumentParser()
for name in ['code', 'tasks', 'results', 'go', 'annotations', 'proteins', 'heldout_train',
             'resources', 'predictor', 'aspect', 'output']:
    ap.add_argument('--' + name.replace('_', '-'), required=True, dest=name)
ap.add_argument('--summary', help='accepted summary JSON to cross-check before/stage1 metrics')
a = ap.parse_args()
sys.path.insert(0, a.code)
from gaef_patches.ontology import Ontology, NAMESPACES, FUNC_DICT
import pipeline.genome_scores_evaluator as gse
from pipeline.genome_scores_evaluator import evaluate_prediction_metrics, load_annotations
from analysis.evaluation_scope import load_temporal_proteins
from analysis.aggregate_revision_audit import checked_reference_labels
from taxon_consistency import adjust_ortools as tc

# ---- capture PR points -------------------------------------------------------
points = []
_orig = gse.evaluate_annotations_optimized


def spy(go, labels, preds, ic_cache=None):
    res = _orig(go, labels, preds, ic_cache)
    points.append((res[1], res[2]))
    return res


gse.evaluate_annotations_optimized = spy


def evaluate(ann, preds):
    points.clear()
    m = evaluate_prediction_metrics(go, a.aspect, ann, preds, reuse_ic=True)
    if m is None:
        return None
    P = np.array([p for p, r in points[-101:]]); R = np.array([r for p, r in points[-101:]])
    idx = np.argsort(R)
    m['aupr_original_rule'] = float(np.trapz(P[idx], R[idx]))
    m['n_eval'] = len(ann)
    return m


# ---- labels, IC: identical to aggregate_revision_audit.aggregate -------------
go = Ontology(a.go)
temporal = load_temporal_proteins(a.proteins)
heldout_train = set(Path(a.heldout_train).read_text().split())
reference_manifest = json.loads((Path(a.annotations) / 'reference-manifest.json').read_text())
reference_records = {r['taxon']: r for r in reference_manifest['records']}
anc_cache = {}


def ancestors(term):
    if term not in anc_cache:
        anc_cache[term] = go.get_ancestors(term)
    return anc_cache[term]


split = 'heldout' if a.predictor.endswith('_heldout') else 'timeset'
counts = Counter()
for path in sorted(Path(a.annotations).glob('annots_taxon_*.tsv')):
    taxon = path.stem.split('_')[-1]
    for protein, terms in load_annotations(path).items():
        closed = set().union(*(ancestors(g) for g in terms))
        if split == 'timeset' and protein not in temporal.get(taxon, ()):
            counts.update(closed)
        if split == 'heldout' and taxon in heldout_train:
            counts.update(closed)
go.calculate_ic_from_counts(counts)

# ---- taxon constraints, as revision_audit_run.run --------------------------
res = Path(a.resources)
in_tax, never_tax = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
in_tax, never_tax = go.normalize_constraints(in_tax), go.normalize_constraints(never_tax)
tax_isa, tax_disjoint, tax_union = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')


def lineage_of(taxon):
    lin = {tc.normalize_taxon_id(taxon)}
    lin.update(tc.get_all_ancestors_of_list(tax_isa, list(lin)))
    for union, members in tax_union.items():
        if members & lin:
            lin.add(union)
    return lin


tasks = [tuple(l.split()) for l in Path(a.tasks).read_text().splitlines() if l.strip()]
tasks = [t for t in tasks if t[0] == a.predictor and t[1] == a.aspect]
conds = ['before', 'stage1', 'stage1_unknown', 'filter_all']
preds = {c: {} for c in conds}
ann = {}
changed = {c: set() for c in conds}
stats = Counter()
removed_terms = Counter()
per_taxon = []
for model, aspect, taxon in tasks:
    data = json.load(gzip.open(Path(a.results) / model / aspect / f'{taxon}.json.gz', 'rt'))
    native = data['native_tau']
    labels = checked_reference_labels(data, a.annotations, reference_records[taxon], temporal)
    by = {r['mode']: r for r in data['runs'] if r['tau'] == native}
    lin = lineage_of(taxon)
    direct = {g for g in set(in_tax) | set(never_tax)
              if any(t not in lin for t in in_tax.get(g, ())) or any(t in lin for t in never_tax.get(g, ()))}
    imp_cache = {}

    def impossible(g):
        if g not in imp_cache:
            imp_cache[g] = bool(ancestors(g) & direct)
        return imp_cache[g]

    before = by['before']['test_scores']
    runs = {'before': before, 'stage1': by['stage1']['test_scores']}
    if 'stage1_unknown' in by:
        runs['stage1_unknown'] = by['stage1_unknown']['test_scores']
    runs['filter_all'] = {p: {g: s for g, s in sc.items() if not impossible(g)} for p, sc in before.items()}
    t = Counter()
    for p in labels:
        b = before.get(p, {})
        lab = set().union(*(ancestors(g) for g in labels[p])) if labels[p] else set()
        f = runs['filter_all'].get(p, {})
        s1 = runs['stage1'].get(p, {})
        for g, s in b.items():
            if g not in f:
                t['filter_removed'] += 1
                t['filter_removed_above_tau' if s > native else 'filter_removed_below_tau'] += 1
                t['filter_removed_true_label'] += int(g in lab)
                removed_terms[g] += 1
            if s1.get(g) != s:
                t['stage1_changed'] += 1
                t['stage1_changed_true_label'] += int(g in lab)
                if g in f:
                    t['stage1_changed_not_in_filter'] += 1
        t['labels_impossible_pairs'] += sum(1 for g in lab if impossible(g))
    for c, sc in runs.items():
        for p in labels:
            key = taxon + ':' + p
            s = sc.get(p, {})
            if c != 'before' and s != before.get(p, {}):
                changed[c].add(key)
            prop = defaultdict(float)
            for term, v in s.items():
                for anc in ancestors(term):
                    prop[anc] = max(prop[anc], v)
            preds[c][key] = dict(prop)
    for p, ts in labels.items():
        ann[taxon + ':' + p] = set(ts)
    t['n_test'] = len(labels)
    t['filter_changed_proteins'] = sum(1 for p in labels if runs['filter_all'].get(p, {}) != before.get(p, {}))
    t['stage1_changed_proteins'] = sum(1 for p in labels if runs['stage1'].get(p, {}) != before.get(p, {}))
    per_taxon.append({'taxon': taxon, **t})
    stats.update(t)
    del data, runs

namespace = go.get_namespace_terms(NAMESPACES[a.aspect]) - {FUNC_DICT[a.aspect]}
valid = {p for p, ts in ann.items() if set().union(*(ancestors(g) for g in ts)) & namespace}
out = {'predictor': a.predictor, 'aspect': a.aspect, 'n_tasks': len(tasks),
       'n_valid_proteins': len(valid), 'stats': stats, 'per_taxon': per_taxon,
       'top_filter_removed_terms': [(g, n, go.get_term(g)['name']) for g, n in removed_terms.most_common(25)],
       'metrics': {}, 'modified_subset': {}}
for c in conds:
    if not preds[c]:
        continue
    out['metrics'][c] = evaluate(ann, preds[c])
    sub = changed[c] & valid
    out['modified_subset'][c] = {'n': len(sub), 'n_taxa': len({p.split(':')[0] for p in sub})}
    if sub and c != 'before':
        out['modified_subset'][c]['before'] = evaluate({p: ann[p] for p in sub}, preds['before'])
        out['modified_subset'][c]['after'] = evaluate({p: ann[p] for p in sub}, preds[c])
    print(c, {k: out['metrics'][c][k] for k in ['fmax', 'smin', 'avg_auc', 'aupr', 'aupr_original_rule']},
          out['modified_subset'][c]['n'], flush=True)
if a.summary and Path(a.summary).exists():
    rows = json.loads(Path(a.summary).read_text())['rows']
    out['crosscheck'] = {r['mode']: {k: r['metrics'][k] for k in ['fmax', 'smin', 'avg_auc', 'aupr']}
                         for r in rows if r['tau'] == (0.3 if a.predictor == 'deepgo-se' else 0.001)}
    print('crosscheck', out['crosscheck'])
Path(a.output).write_text(json.dumps(out, indent=1, default=float))
print(dict(stats))
