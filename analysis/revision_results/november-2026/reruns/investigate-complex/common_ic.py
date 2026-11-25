"""Shared loaders for the Stage 2 investigation (mirrors revision_audit_run)."""
import gzip, json, sys
from collections import defaultdict, Counter
from pathlib import Path

TAU = {'deepgo-se': .3, 'seq_sim_timeset': .001, 'sprof-go': .001}
PREDS = ['seq_sim_timeset', 'deepgo-se', 'sprof-go']


def setup(code):
    sys.path.insert(0, code)


def resources(code, go_file):
    setup(code)
    from gaef_patches.ontology import Ontology
    from complex_coherence import adjust_ortools as cc
    go = Ontology(go_file)
    hierarchy = defaultdict(set, {g: set(o['is_a']) for g, o in go.ont.items()})
    isa = Ontology(go_file, with_rels=False)
    isa_h = {g: set(o['is_a']) for g, o in isa.ont.items()}
    homomers = go.normalize_terms(cc.load_homodimer_terms(Path(code)/'complex_coherence/protein_complexes.tsv'))
    complexes = cc.get_heteromeric_complexes(homomers, set(), isa_h)
    return go, hierarchy, isa_h, complexes, cc


def stage_file(accepted, pred, mode, taxon):
    fs = sorted((Path(accepted)/pred/'cc'/mode).glob(f'optimized_fold_*_taxon_{taxon}.tsv'))
    if len(fs) != 1:
        raise FileNotFoundError(f'{pred} {mode} {taxon}: {fs}')
    return fs[0]


def active(scores, tau):
    return {p: {g for g, s in sc.items() if s > tau} for p, sc in scores.items()}


def load_cp_index(path, go, complexes):
    idx = json.load(open(path))
    out = defaultdict(set)
    for key, accs in idx.items():
        taxon, term = key.split('::', 1)
        term = go.canonical_id(term)
        if term in complexes:
            out[taxon+'::'+term].update(accs)
    return out


def load_uniprot_index(path):
    with gzip.open(path, 'rt') as fh:
        return json.load(fh)


def protein_to_accession(pid, mapping):
    from analysis.complex_validation.validate_stage2 import protein_to_accession as f
    return f(pid, mapping)
