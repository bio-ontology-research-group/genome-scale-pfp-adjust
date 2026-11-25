"""Diagnose accepted Stage 2 repairs per singleton complex (one predictor x organism).

For every organism-complex singleton in the Stage 1 output: the lone protein
and score, the H3 candidate pool, what standard/demotion-only Stage 2 did,
local demotion and promotion costs (the solver's integer costs), and where
Complex Portal subunits and UniProt-annotated members rank by predicted score.
"""
import argparse, gzip, json
from collections import defaultdict
from pathlib import Path
import common_ic as C


def load_active_complex(path, tau, complexes):
    out = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            for tok in cols[1:]:
                g, s = tok.split('|')
                if g in complexes and float(s) > tau:
                    out[cols[0]].add(g)
    return out


def icost(c):
    return max(1, int(c * 1e6)) / 1e6


def main():
    ap = argparse.ArgumentParser()
    for k in ['code', 'go', 'accepted', 'cp_index', 'uniprot_dir', 'idmap', 'predictor', 'taxon', 'out']:
        ap.add_argument('--'+k, required=True)
    a = ap.parse_args()
    go, H, isa_h, complexes, cc = C.resources(a.code, a.go)
    tau = C.TAU[a.predictor]
    scores = cc.load_predictions(C.stage_file(a.accepted, a.predictor, 'stage1', a.taxon))
    ann = C.active(scores, tau)
    std = load_active_complex(C.stage_file(a.accepted, a.predictor, 'stage2_standard', a.taxon), tau, complexes)
    dem = load_active_complex(C.stage_file(a.accepted, a.predictor, 'stage2_demotion', a.taxon), tau, complexes)
    members = defaultdict(set)
    for p, ts in ann.items():
        for g in ts & complexes:
            members[g].add(p)
    singles = sorted(g for g, m in members.items() if len(m) == 1)
    cp = C.load_cp_index(a.cp_index, go, complexes)
    mapping = json.load(open(a.idmap))
    proteome = set(scores)
    inv = defaultdict(set)
    for p in proteome:
        acc = C.protein_to_accession(p, mapping)
        if acc:
            inv[acc].add(p)
    upath = Path(a.uniprot_dir)/f'{a.taxon}.json.gz'
    uni = C.load_uniprot_index(upath)['index'] if upath.exists() else {}
    anc_cache, desc_cache = {}, {}
    def anc(g):
        if g not in anc_cache: anc_cache[g] = cc.get_ancestors(g, H)
        return anc_cache[g]
    def desc(g):
        if g not in desc_cache: desc_cache[g] = cc.get_descendants(g, H)
        return desc_cache[g]
    def prom_cost(q, X):
        return sum(icost(abs(scores[q].get(g, 0.0) - tau)) for g in {X} | anc(X) if g not in ann[q] and g not in cc.ROOT_TERMS)
    def dem_cost(L, X):
        return sum(icost(abs(scores[L].get(g, 0.0) - tau)) for g in {X} | desc(X) if g in ann[L])
    recs = []
    names = sorted(proteome)
    for X in singles:
        L = next(iter(members[X]))
        order = sorted(names, key=lambda p: (-scores[p].get(X, 0.0), p))
        rank = {p: i+1 for i, p in enumerate(order)}
        top50 = set(order[:50])
        others = [p for p in order if p != L]
        pos = [p for p in others if scores[p].get(X, 0.0) > 0]
        cand = [p for p in order[:50] if p != L]
        pc = {q: prom_cost(q, X) for q in cand}
        best = min(cand, key=lambda q: (pc[q], q)) if cand else None
        promoted = sorted(p for p, ts in std.items() if X in ts and p != L)
        r = dict(term=X, lone=L, lone_score=scores[L].get(X, 0.0), n_proteins=len(names),
                 n_pos_inactive=len(pos), max_other_score=scores[others[0]].get(X, 0.0) if others else None,
                 n_top50_pos=sum(scores[q].get(X, 0.0) > 0 for q in cand),
                 n_anc=len(anc(X) - cc.ROOT_TERMS), n_active_desc_on_lone=len((desc(X) & ann[L])),
                 dem_cost=dem_cost(L, X), best_prom_cost=pc.get(best), best_prom=best,
                 n_cand_at_min_cost=sum(pc[q] == pc[best] for q in cand) if cand else 0,
                 std_lone_kept=X in std.get(L, set()), std_promoted=promoted,
                 std_promoted_rank=[rank[p] for p in promoted], std_promoted_score=[scores[p].get(X, 0.0) for p in promoted],
                 std_promoted_cost=[prom_cost(p, X) for p in promoted],
                 dem_lone_kept=X in dem.get(L, set()))
        key = f'{a.taxon}::{X}'
        if key in cp:
            subs = []
            for acc in sorted(cp[key]):
                for p in sorted(inv.get(acc, ())):
                    subs.append(dict(acc=acc, entry=p, score=scores[p].get(X, 0.0), rank=rank[p], in_top50=p in top50,
                                     is_lone=p == L, prom_cost=None if p == L else prom_cost(p, X)))
            r['cp'] = dict(n_listed=len(cp[key]), n_mapped=len({s['acc'] for s in subs}), subunits=subs,
                           lone_listed=C.protein_to_accession(L, mapping) in cp[key])
        um = sorted(p for p in names if X in uni.get(p, {}))
        r['uni'] = dict(lone_class=uni.get(L, {}).get(X),
                        members=[dict(entry=p, cls=uni[p][X], score=scores[p].get(X, 0.0), rank=rank[p], in_top50=p in top50,
                                      prom_cost=prom_cost(p, X)) for p in um if p != L])
        recs.append(r)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(a.out, 'wt') as fh:
        json.dump(dict(predictor=a.predictor, taxon=a.taxon, tau=tau, n_singletons=len(singles),
                       uniprot_index=upath.exists(), records=recs), fh)
    print(a.predictor, a.taxon, len(singles))


if __name__ == '__main__':
    main()
