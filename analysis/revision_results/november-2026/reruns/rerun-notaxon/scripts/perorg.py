"""Per-organism summary of rerun outputs (no test scores): demotions, true-lineage consistency, fixed taxa."""
import gzip, json, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
root = Path(sys.argv[1]); tasks = [l.split() for l in Path(sys.argv[2]).read_text().splitlines() if l.strip()]
def one(t):
    pred, asp, tax = t
    d = json.load(gzip.open(root / pred / asp / f'{tax}.json.gz', 'rt'))
    rec = {'predictor': pred, 'aspect': asp, 'taxon': tax, 'runs': {}}
    for r in d['runs']:
        rec['runs'][r['mode']] = {k: r.get(k) for k in ['demotions', 'known_taxon_consistent', 'gaef_taxon_consistent', 'fixed_taxa',
                                                        'lineage_inference', 'solver_taxons', 'seconds_solve', 'active_predictions']}
        rec['runs'][r['mode']]['n_violating'] = len(r['violating_terms'])
        rec['runs'][r['mode']]['n_modified_test'] = len(r.get('modified_test_proteins', []))
    return rec
with ProcessPoolExecutor(int(sys.argv[4]) if len(sys.argv) > 4 else 4) as ex:
    out = list(ex.map(one, tasks))
Path(sys.argv[3]).write_text(json.dumps(out))
print(len(out))
