"""Rebuild local baseline labels while preserving recorded protein/taxon splits."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gaef_patches.ontology import Ontology, NAMESPACES


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def prepare(args):
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Use a new baseline-input directory')
    refs = Path(args.references)
    manifest_path = refs/'reference-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest['go_sha256'] != digest(args.go):
        raise ValueError('Reference ontology differs')
    records = {r['taxon']: r for r in manifest['records']}
    split_rows = {}
    needed = set()
    for split, source in [('train', args.train), ('val', args.validation)]:
        with Path(source).open() as handle:
            split_rows[split] = list(csv.reader(handle, delimiter='\t'))
        if any(len(row) < 2 for row in split_rows[split]):
            raise ValueError('Malformed training split')
        needed.update(row[0] for row in split_rows[split])
    labels = {}
    for taxon, row in records.items():
        path = refs/row['reference_file']
        if digest(path) != row['reference_sha256']:
            raise ValueError(f'Reference checksum differs: {path}')
        labels[taxon] = {}
        with path.open() as handle:
            for line in handle:
                protein, *terms = line.rstrip('\n').split('\t')
                if protein in needed:
                    labels[taxon][protein] = {t for t in terms if t}

    go = Ontology(args.go)
    output.mkdir(parents=True, exist_ok=True)
    source_hashes = {str(manifest_path):digest(manifest_path), str(args.go):digest(args.go)}
    files = {}
    split_counts = {}
    seen = set()
    test_taxa = set(Path(args.test_taxa).read_text().split())
    source_hashes[str(args.test_taxa)] = digest(args.test_taxa)
    for split, source in [('train',args.train),('val',args.validation)]:
        source_hashes[str(source)] = digest(source)
        target = output/(split+'_proteins.tsv')
        counts = dict(original=0, retained=0, empty_removed=0, changed=0)
        with target.open('w') as dest:
            for row in split_rows[split]:
                if len(row) < 2:
                    raise ValueError('Malformed training split')
                protein, taxon, *old = row
                if protein in seen or taxon in test_taxa:
                    raise ValueError(f'Duplicate/train-validation/test overlap: {protein}/{taxon}')
                seen.add(protein)
                if taxon not in labels or protein not in labels[taxon]:
                    raise ValueError(f'Training protein missing from source: {protein}/{taxon}')
                terms = set(labels[taxon][protein])
                counts['original'] += 1
                counts['changed'] += go.normalize_terms(old) != terms
                if not terms:
                    counts['empty_removed'] += 1
                    continue
                dest.write('\t'.join([protein,taxon]+sorted(terms))+'\n')
                counts['retained'] += 1
        if not counts['retained']:
            raise ValueError('Empty corrected training split')
        split_counts[split] = counts
        files[target.name] = digest(target)
    vocabularies = {}
    for aspect, namespace in NAMESPACES.items():
        source = Path(args.terms_root)/aspect/'terms.pkl'
        source_hashes[str(source)] = digest(source)
        old = pd.read_pickle(source)['gos'].values.flatten()
        terms = set().union(*(go.get_ancestors(t) for t in old)) & go.get_namespace_terms(namespace)
        if not terms:
            raise ValueError(f'Empty vocabulary: {aspect}')
        folder = output/aspect
        folder.mkdir()
        target = folder/'terms.pkl'
        pd.DataFrame({'gos':sorted(terms)}).to_pickle(target)
        files[str(target.relative_to(output))] = digest(target)
        vocabularies[aspect] = len(terms)
    result = dict(schema_version=1, definition='Source-rebuilt labels; original split membership; remove empty-label proteins; canonical inherited vocabulary with valid ancestors',
                  sources=source_hashes, files=files, split_counts=split_counts,
                  vocabularies=vocabularies, go_sha256=digest(args.go),
                  code_sha256=digest(__file__))
    (output/'training-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['references','go','train','validation','test-taxa','terms-root','output']:
        parser.add_argument('--'+name, required=True)
    print(json.dumps(prepare(parser.parse_args()),indent=2))
