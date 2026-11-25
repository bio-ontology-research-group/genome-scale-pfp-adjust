"""Per-sub-ontology lineage outcomes under the repair cost (Data 2) and the lineage objective (lineage_variants.py records).
Writes a LaTeX tabular body and a per-organism TSV."""
import csv, json, glob, sys, collections, pathlib
lv, data2, out_tex, out_tsv = sys.argv[1:5]
RANKS = ['domain','phylum','class','order','family','genus','species']
SETS = [('deepgo-se','DGO-E'),('seq_sim_timeset','SS-E'),('sprof-go','SP-E'),('seq_sim_heldout','SS-H'),('mlp_heldout','MLP-H')]
ASP = [('cc','CC'),('mf','MF'),('bp','BP'),('joint','Joint')]
cat = collections.defaultdict(collections.Counter); deep = collections.defaultdict(collections.Counter)
rows = []
for r in csv.DictReader(open(data2), delimiter='\t'):
    n = int(r['supported_n_inferred_nodes'])
    c = 'unresolved' if n == 0 else ('correct' if r['supported_is_exact_prefix'] == 'True' else 'wrong')
    k = ('repair cost', r['predictor'], r['aspect'])
    cat[k][c] += 1; deep[k][r['supported_deepest_correct_rank']] += 1
    rows.append(['repair cost', r['predictor'], r['taxon'], r['aspect'], c, r['supported_deepest_correct_rank'], r['supported_first_divergence_rank'], ''])
for sub in ['lv','lv-heldout']:
    for f in sorted(glob.glob(f'{lv}/{sub}/*.json')):
        pred, taxon = pathlib.Path(f).stem.rsplit('_',1)
        for r in json.load(open(f))['records']:
            if r['variant'] != 'siglogit': continue
            assert r['status'] == 'OPTIMAL' and not r.get('n_undetermined')
            k = ('lineage objective', pred, r['aspect'])
            cat[k][r['category']] += 1; deep[k][r['deepest']] += 1
            inf = ';'.join(f"{rk}={','.join(sorted(r['inferred_ranks'][rk]))}" for rk in RANKS if r['inferred_ranks'].get(rk))
            rows.append(['lineage objective', pred, taxon, r['aspect'], r['category'], r['deepest'], r['divergence'], inf])
lines = []
for p, lab in SETS:
    for a, alab in ASP:
        cells = []
        for obj in ['repair cost','lineage objective']:
            c, d = cat[(obj,p,a)], deep[(obj,p,a)]
            n = sum(c.values()); assert n in (50, 219), (obj,p,a,n)
            assert not any(d.get(x) for x in ['family','genus','species'])
            cells += [f"{c['correct']}/{c['wrong']}/{c['unresolved']}", f"{d['domain']}/{d['phylum']}/{d['class']}/{d['order']}"]
        lines.append(f"{lab} & {alab} & " + ' & '.join(cells) + r' \\')
    if lab != 'MLP-H': lines.append(r'\addlinespace')
pathlib.Path(out_tex).write_text('\n'.join(lines) + '\n')
with open(out_tsv, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(['objective','predictor','taxon','aspect','outcome','deepest_correct_rank','first_divergence_rank','inferred_taxa'])
    w.writerows(sorted(rows))
print('\n'.join(lines)); print(len(rows), 'rows')
