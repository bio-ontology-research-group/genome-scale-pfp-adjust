"""Build verified UniProt entry-name/accession mappings from idmapping files."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def build(paths):
    mapping = {}
    for path in paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as fh:
            for line in fh:
                fields = line.rstrip('\n').split('\t')
                if len(fields) == 3 and fields[1] == 'UniProtKB-ID':
                    accession, _, entry = fields
                    if entry in mapping and mapping[entry] != accession:
                        raise ValueError(f'Ambiguous entry name: {entry}')
                    mapping[entry] = accession
    return mapping


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    mapping = build(args.inputs)
    if not mapping:
        raise ValueError('No UniProtKB-ID mappings found')
    Path(args.output).write_text(json.dumps(mapping, sort_keys=True) + '\n')
