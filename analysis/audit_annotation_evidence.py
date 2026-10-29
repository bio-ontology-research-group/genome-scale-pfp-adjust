"""Audit missing complex partners by evidence set, using complete proteome DATs.

This is an annotation-coverage sensitivity analysis, not an accuracy benchmark:
computational annotations are not treated as independent biological truth.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaef_patches.ontology import Ontology
from complex_coherence.adjust_ortools import load_homodimer_terms, get_heteromeric_complexes

DIRECT = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
CURATED = DIRECT | {'TAS', 'IC'}


def file_sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def entries(path):
    entry, accession, annotations = None, None, []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('ID   '):
                entry = line.split()[1]
            elif line.startswith('AC   ') and accession is None:
                accession = line[5:].split(';')[0].strip()
            elif line.startswith('DR   GO;'):
                fields = [x.strip() for x in line[5:].split(';')]
                if len(fields) >= 4:
                    evidence = fields[3].split(':')[0]
                    annotations.append((fields[1], evidence))
            elif line.startswith('//'):
                if not entry or not accession:
                    raise ValueError(f'Missing entry/accession in {path}')
                yield entry, accession, annotations
                entry, accession, annotations = None, None, []
    if entry is not None:
        raise ValueError(f'Unterminated entry in {path}')


def summarize(counts):
    present = {g for g, n in counts.items() if n}
    singleton = {g for g in present if counts[g] == 1}
    return {'present': len(present), 'singletons': len(singleton),
            'coherence': 1-len(singleton)/len(present) if present else None,
            'singleton_terms': sorted(singleton)}


def run(args):
    go = Ontology(args.go)
    isa_go = Ontology(args.go, with_rels=False)
    hierarchy = {g: set(obj['is_a']) for g, obj in isa_go.ont.items()}
    complexes = get_heteromeric_complexes(go.normalize_terms(load_homodimer_terms(ROOT/'complex_coherence/protein_complexes.tsv')), set(), hierarchy)
    taxa = set(Path(args.taxa).read_text().split())
    with open(args.proteomes) as fh:
        index = {r['TaxonID']: r for r in csv.DictReader(fh, delimiter='\t')}
    mapping, rows = {}, []
    cache = {}
    for taxon in sorted(taxa):
        row = index[taxon]
        paths = sorted(Path(args.root).glob(f"*/{row['ProteomeID']}/{row['ProteomeID']}_{taxon}.dat.gz"))
        if len(paths) != 1:
            raise ValueError(f'Expected exactly one DAT for {taxon}: {paths}')
        counts = {key: Counter() for key in ['direct', 'curated', 'all_same_proteins', 'all_proteome']}
        n_proteins = Counter()
        evidence_counts = Counter()
        for entry, acc, annotations in entries(paths[0]):
            if entry in mapping and mapping[entry] != acc:
                raise ValueError(f'Conflicting entry name {entry}')
            mapping[entry] = acc
            terms = {key: set() for key in ['direct', 'curated', 'all']}
            for term, evidence in annotations:
                evidence_counts[evidence] += 1
                if term not in cache:
                    cache[term] = go.get_ancestors(term) & complexes
                terms['all'].update(cache[term])
                if evidence in DIRECT:
                    terms['direct'].update(cache[term])
                if evidence in CURATED:
                    terms['curated'].update(cache[term])
            has_curated = any(e in CURATED for _, e in annotations)
            for key in ['direct', 'curated']:
                counts[key].update(terms[key])
                n_proteins[key] += bool(terms[key])
            counts['all_proteome'].update(terms['all'])
            if has_curated:
                counts['all_same_proteins'].update(terms['all'])
            n_proteins['all_proteome'] += 1
            n_proteins['with_curated_any_aspect'] += has_curated
        summaries = {key: summarize(value) for key, value in counts.items()}
        initial = set(summaries['curated']['singleton_terms'])
        rows.append({'taxon': taxon, 'file': str(paths[0]),
                     'sha256': file_sha(paths[0]),
                     'protein_counts': dict(n_proteins), 'evidence_counts': dict(evidence_counts),
                     'complex_counts': summaries,
                     'curated_singletons_with_partner_all_same_proteins': sum(counts['all_same_proteins'][g] >= 2 for g in initial),
                     'curated_singletons_with_partner_all_proteome': sum(counts['all_proteome'][g] >= 2 for g in initial)})
        print(taxon, n_proteins['all_proteome'], summaries['curated']['singletons'], flush=True)
    Path(args.mapping).write_text(json.dumps(mapping, sort_keys=True)+'\n')
    Path(args.output).write_text(json.dumps({'provenance': {'script_sha256': file_sha(__file__), 'go_sha256': file_sha(args.go), 'complex_exclusions_sha256': file_sha(ROOT/'complex_coherence/protein_complexes.tsv'), 'proteome_index_sha256': file_sha(args.proteomes)}, 'direct_evidence': sorted(DIRECT), 'curated_evidence': sorted(CURATED),
        'n_complex_terms': len(complexes), 'rows': rows}, indent=2)+'\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'proteomes', 'taxa', 'go', 'mapping', 'output']:
        ap.add_argument('--'+name, required=True)
    run(ap.parse_args())
