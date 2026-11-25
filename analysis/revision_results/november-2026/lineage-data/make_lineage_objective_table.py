"""Supplementary Table S13 for the lineage objective (siglogit records of lineage_variants.py)."""
import json, glob, sys, collections, pathlib
lv, out = sys.argv[1:3]
RANKS = ['none','domain','phylum','class','order','family','genus','species']
R7 = RANKS[1:]
SETS = [('deepgo-se','DGO-E'),('seq_sim_timeset','SS-E'),('sprof-go','SP-E'),('seq_sim_heldout','SS-H'),('mlp_heldout','MLP-H')]
ASP = [('cc','CC'),('mf','MF'),('bp','BP'),('joint','JOINT')]
deep = collections.defaultdict(collections.Counter); div = collections.defaultdict(collections.Counter)
per = collections.defaultdict(dict)
for sub in ['lv','lv-heldout']:
    for f in sorted(glob.glob(f'{lv}/{sub}/*.json')):
        pred, taxon = pathlib.Path(f).stem.rsplit('_',1)
        for r in json.load(open(f))['records']:
            if r['variant'] != 'siglogit': continue
            assert r['status'] == 'OPTIMAL' and not r.get('n_undetermined')
            deep[(pred,r['aspect'])][r['deepest']] += 1; div[(pred,r['aspect'])][r['divergence']] += 1
            per[(pred,taxon)][r['aspect']] = r['inferred_ranks']
L = [r'\begingroup\scriptsize', r'\setlength{\tabcolsep}{4pt}', r'\begin{longtable}{lrrrrrrrrr}',
     r'\caption{Lineage inference with the lineage objective: number of organisms by deepest correct rank (Deepest) and by first divergent rank (Divergence). A taxon counts as inferred only if every optimal solution requires it. None denotes no resolved correct rank or no observed divergence, respectively. SS: Seq-Sim; DGO: DeepGO-SE; SP: SPROF-GO; MLP: MLP-ESM2. E: entry-date split ($n=50$); H: heldout split ($n=219$).}\label{tab:supp-lineage}\\']
head = [r'\toprule', r'Model/aspect & Measure & None & Domain & Phylum & Class & Order & Family & Genus & Species \\', r'\midrule']
L += head + [r'\endfirsthead'] + head + [r'\endhead']
for p, lab in SETS:
    for a, alab in ASP:
        for name, c in [('Deepest', deep[(p,a)]), ('Divergence', div[(p,a)])]:
            assert sum(c.values()) in (50, 219) and set(c) <= set(RANKS)
            L.append(f"{lab} {alab} & {name} & " + ' & '.join(str(c.get(r,0)) for r in RANKS) + r' \\')
L += [r'\bottomrule', r'\end{longtable}', r'\endgroup']
agree = []; rows = []
for p, lab in SETS:
    orgs = [k for k in per if k[0] == p]; ok = 0; cells = []
    for o in orgs:
        good = True
        for rk in R7:
            sets = {frozenset(per[o][a].get(rk, [])) for a in ['cc','mf','bp']} - {frozenset()}
            if len(sets) > 1: good = False
        ok += good
    for rk in R7:
        comp = ident = 0
        for o in orgs:
            raw = [frozenset(per[o][a].get(rk, [])) for a in ['cc','mf','bp']]
            ne = {s for s in raw if s}
            comp += len(ne) <= 1; ident += (all(raw) and len(ne) == 1)
        cells.append(f"{comp}/{ident}")
    rows.append(f"{lab} & " + ' & '.join(cells) + f" & {ok}/{len(orgs)}" + r" \\")
    agree.append(f"{lab} {ok}/{len(orgs)}")
L.append(r'\par Agreement of CC, MF and BP by rank: organisms without conflicting inferred taxa / organisms in which all three sub-ontologies infer the same taxon. All ranks: organisms without a conflict at any rank.')
L += [r'\begin{center}\scriptsize\begin{tabular}{lrrrrrrrr}\toprule', r'Model & Domain & Phylum & Class & Order & Family & Genus & Species & All ranks \\', r'\midrule'] + rows + [r'\bottomrule\end{tabular}\end{center}']
pathlib.Path(out).write_text('\n'.join(L) + '\n'); print('\n'.join(L[-48:]))
