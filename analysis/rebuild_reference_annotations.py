"""Rebuild reference labels from source evidence using only is_a/part_of closure."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from analysis.audit_annotation_evidence import CURATED, entries
from gaef_patches.ontology import Ontology
from pipeline.genome_scores_evaluator import load_annotations


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as fh:
        for block in iter(lambda: fh.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def rebuild(args):
    out = Path(args.output)
    if (out/'reference-manifest.json').exists():
        raise FileExistsError('Use a new output directory for a new reference build')
    originals = sorted(Path(args.original_annotations).glob('annots_taxon_*.tsv'))
    if not originals:
        raise ValueError('No original annotation files identify the reference scope')
    with Path(args.proteome_index).open() as handle:
        index = {r['TaxonID']: r for r in csv.DictReader(handle, delimiter='\t')}
    sources = []
    for original in originals:
        taxon = original.stem.removeprefix('annots_taxon_')
        proteome = index[taxon]['ProteomeID']
        matches = list(Path(args.source_root).glob(f'*/{proteome}/{proteome}_{taxon}.dat.gz'))
        if len(matches) != 1:
            raise ValueError(f'Expected one source DAT for {taxon}: {matches}')
        sources.append((taxon, original, matches[0]))
    go = Ontology(args.go)
    cache = {}
    def closure(term):
        if term not in cache:
            cache[term] = go.get_ancestors(term)
        return cache[term]
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for taxon, original, source in sources:
        previous = load_annotations(original)
        seen = set()
        counts = dict(proteins=0, positive_reference_proteins=0, reference_pairs=0,
                      old_extra_pairs=0, new_pairs_absent_from_old=0,
                      unknown_or_obsolete_source_annotations=0)
        uncovered = []
        target = out/original.name
        partial = target.with_suffix('.tsv.part')
        with partial.open('w') as fh:
            for protein, accession, annotations in entries(source):
                if protein in seen:
                    raise ValueError(f'Duplicate source entry {protein} in {source}')
                seen.add(protein)
                direct = {g for g, evidence in annotations if evidence in CURATED}
                labels = set().union(*(closure(g) for g in direct))
                old_labels = set().union(*(closure(g) for g in previous.get(protein, ())))
                counts['proteins'] += 1
                counts['positive_reference_proteins'] += bool(labels)
                counts['reference_pairs'] += len(labels)
                counts['old_extra_pairs'] += len(old_labels-labels)
                counts['new_pairs_absent_from_old'] += len(labels-old_labels)
                counts['unknown_or_obsolete_source_annotations'] += sum(go.canonical_id(g) is None for g in direct)
                if labels and protein not in previous:
                    uncovered.append(protein)
                fh.write(protein+'\t'+'\t'.join(sorted(labels))+'\n')
        partial.replace(target)
        records.append(dict(taxon=taxon, source_file=str(source), source_sha256=digest(source),
                            original_file=str(original), original_sha256=digest(original),
                            reference_file=target.name, reference_sha256=digest(target),
                            new_reference_proteins_not_in_original_rows=uncovered,
                            original_reference_proteins_not_in_source=sorted(set(previous)-seen),
                            **counts))
        print(taxon, counts['positive_reference_proteins'], counts['old_extra_pairs'], len(uncovered), flush=True)
    manifest = dict(schema_version=1, definition='positive reference evidence; canonical primary GO IDs; is_a and part_of closure only',
                    evidence_codes=sorted(CURATED), go_sha256=digest(args.go),
                    proteome_index_sha256=digest(args.proteome_index),
                    code_sha256={str(p.relative_to(ROOT)):digest(p) for p in
                                 [Path(__file__).resolve(), ROOT/'analysis/audit_annotation_evidence.py', ROOT/'gaef_patches/ontology.py']},
                    records=records)
    (out/'reference-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source-root', 'proteome-index', 'original-annotations', 'go', 'output']:
        parser.add_argument('--'+name, required=True)
    rebuild(parser.parse_args())
