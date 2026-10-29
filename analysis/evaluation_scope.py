"""Explicit protein membership for the entry-date benchmark."""
import csv
import hashlib
from pathlib import Path


def load_temporal_proteins(path):
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            for taxon in row['taxon_id'].split(';'):
                result.setdefault(taxon, set()).add(row['protein_id'])
    return result


def select_annotations(annotations, taxon, temporal_proteins):
    """Use exact UniProt entry names; do not infer accessions by string slicing."""
    if temporal_proteins is None:
        return annotations
    allowed = temporal_proteins.get(str(taxon), set())
    return {p: terms for p, terms in annotations.items() if p in allowed}


def scope_manifest(path):
    path = Path(path)
    members = load_temporal_proteins(path)
    return {'definition': 'Swiss-Prot entry integration date after 2024-05-23',
            'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'n_proteins': sum(map(len, members.values())), 'n_taxa': len(members)}
