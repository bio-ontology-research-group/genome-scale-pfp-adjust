"""Download the published FunTaxIS-lite server's precomputed constraints.

No constraints are inferred locally. Cache the exact responses and checksum
each table so the baseline can be rerun without querying a changing service.
"""
import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, urljoin
from urllib.request import urlopen

BASE = 'https://www.medcomp.medicina.unipd.it/funtaxis-lite/'


def download(taxa_file, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    manifest = {'source': BASE, 'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'taxa': {}}
    for taxon in sorted(set(Path(taxa_file).read_text().split())):
        query = output / f'{taxon}.response.json'
        if not query.exists():
            with urlopen(BASE + 'searchdata', data=urlencode({'search_word': taxon}).encode(), timeout=60) as r:
                query.write_bytes(r.read())
            time.sleep(.2)
        html = json.loads(query.read_text())['data']
        links = re.findall(r'href="([^"]*/static/constraints_lists/[^"]+)"', html)
        if len(links) != 1:
            raise ValueError(f'Expected one published constraint list for taxon {taxon}; found {links}')
        url = urljoin(BASE, links[0])
        table = output / Path(links[0]).name
        if not table.exists():
            with urlopen(url, timeout=60) as r:
                table.write_bytes(r.read())
        terms = sorted(set(re.findall(r'GO:\d{7}', table.read_text())))
        if not terms:
            raise ValueError(f'No GO terms in {table}')
        manifest['taxa'][taxon] = {
            'url': url, 'file': table.name, 'sha256': hashlib.sha256(table.read_bytes()).hexdigest(),
            'n_terms': len(terms), 'terms': terms,
        }
        (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
        print(taxon, len(terms), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--taxa', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    download(args.taxa, args.output)
