"""Diagnose 'wrong' contig lineages: required taxa for whole genome and sampled 1 Mb contigs,
with and without the true taxon fixed."""
import csv, random, sys, argparse
from collections import defaultdict
from pathlib import Path
ap = argparse.ArgumentParser()
for k in ['code', 'resources', 'predictions', 'cds', 'taxon', 'aspect']:
    ap.add_argument('--' + k, required=True)
ap.add_argument('--genome-size', type=int, required=True)
ap.add_argument('--tau', type=float, default=.001)
a = ap.parse_args()
sys.path.insert(0, a.code)
from gaef_patches.ontology import Ontology, NAMESPACES
from taxon_consistency import adjust_ortools as tc
from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage, RANKS_7
res = Path(a.resources)
go = Ontology(res / 'go-basic.obo')
hier = defaultdict(set, {g: set(o['is_a']) for g, o in go.ont.items()})
raw = go.propagate_scores(go.normalize_scores(tc.load_predictions(a.predictions), NAMESPACES[a.aspect], {}), {})
it, nt = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
it, nt = go.normalize_constraints(it), go.normalize_constraints(nt)
ti, td, tu = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
pm, rm, nm = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
ranked = {t for t, r in rm.items() if r in RANKS_7}
true = get_true_lineage(a.taxon, pm)
print('true lineage sample', list(true)[:5])
def name(t):
    try: return nm.get(t) or nm.get(str(t).replace('NCBITaxon_', ''))
    except Exception: return ''
cds = [(int(r[0]), int(r[1]), r[3]) for r in list(csv.reader(open(a.cds), delimiter='\t'))[1:]]
def run(prots, label):
    c = {p: raw[p] for p in prots}
    ann = {p: {g for g, s in sc.items() if s > a.tau} for p, sc in c.items()}
    costs, _ = tc.compute_demotion_flip_cost(c, a.tau)
    for fix in [None, a.taxon]:
        unc = {}
        adj, taxa, flips = tc.solve_sat_ortools_taxon_assignment(ann, it, nt, hier, costs, ti, td, tu,
            genome_taxon_id=fix, uncertainty=unc, query_taxa=ranked, time_limit=120)
        dem = [(p, g) for p in ann for g in ann[p] if not adj.get(p, {}).get(g, False)]
        cost = sum(costs.get(g, 0) if not isinstance(costs.get(g), dict) else 0 for _, g in dem)
        req = unc.get('required', [])
        print(label, 'fixed=' + str(fix), 'status', unc.get('status'), 'n_prot', len(ann), 'demoted', len(dem),
              'required', [(t, rm.get(t), name(t)) for t in req][:8], flush=True)
        from collections import Counter
        print('   top demoted terms', Counter(g for _, g in dem).most_common(6))
run(sorted(p for p in raw if any(n == p for _, _, n in cds)), 'whole')
rng = random.Random(42)
for L in [1000, 10000, 100000]:
    for i in range(1000): rng.randint(1, max(1, a.genome_size - L))
for i in range(2):
    s = rng.randint(1, a.genome_size - 1000000); e = s + 1000000
    run(sorted({n for x, y, n in cds if x <= e and y >= s and n in raw}), f'1Mb#{i}')
