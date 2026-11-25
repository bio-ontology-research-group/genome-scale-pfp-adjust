"""Recompute the joint (CC+MF+BP) siglogit lineage solve and keep the SELECTED assignment.

lineage_variants.py saves only the certified (required) taxa; the witness
assignment of the first CP-SAT solve (num_search_workers=1, random_seed=42,
sorted variables) is discarded. This recomputes exactly that first solve with
the same inputs (termstats summaries, siglogit costs, joint aspect) and records
the selected non-union lineage chain, its deepest node, and how it relates to
the certified prefix and to the true lineage. The uncertainty checks run after
the first solve and do not change its witness, so they are skipped here; the
saved record's demoted term count and cost are checked for reproduction.
"""
import argparse
import gzip
import json
import sys
import time
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    for a in ['code', 'resources', 'termstats', 'predictor', 'taxa', 'lv-dir', 'lv-scripts', 'output']:
        ap.add_argument('--' + a, required=True)
    ap.add_argument('--variant', default='siglogit')
    ap.add_argument('--time-limit', type=float, default=120)
    a = ap.parse_args()
    sys.path.insert(0, a.code)
    sys.path.insert(0, a.lv_scripts)
    from lineage_variants import variant_inputs
    from gaef_patches.ontology import Ontology
    from taxon_consistency import adjust_ortools as tc
    from analysis.evaluate_lineage_inference import parse_ncbitaxon_obo, get_true_lineage, RANKS_7
    res = Path(a.resources)
    go = Ontology(res / 'go-basic.obo')
    hierarchy = {g: set(o['is_a']) for g, o in go.ont.items()}
    only, never = tc.load_constraints(res / 'go_taxon_constraints_extracted_obo.tsv')
    only, never = go.normalize_constraints(only), go.normalize_constraints(never)
    isa, dis, uni = tc.load_taxon_hierarchy(res / 'taxon_hierarchy.tsv', res / 'ncbitaxon_hierarchy.tsv')
    parents, ranks, names = parse_ncbitaxon_obo(str(res / 'ncbitaxon.obo'))
    constrained = {g: (frozenset(only.get(g, ())), frozenset(never.get(g, ()))) for g in set(only) | set(never)}
    out = Path(a.output)
    done = json.loads(out.read_text()) if out.exists() else {}
    for taxon in a.taxa.split(','):
        if taxon in done:
            continue
        lv = json.loads((Path(a.lv_dir) / f'{a.predictor}_{taxon}.json').read_text())
        saved = [r for r in lv['records'] if r['variant'] == a.variant and r['aspect'] == 'joint']
        assert len(saved) == 1, (taxon, len(saved))
        saved = saved[0]
        stats = {}
        for asp in ['cc', 'mf', 'bp']:
            with gzip.open(Path(a.termstats) / a.predictor / asp / f'{taxon}.json.gz', 'rt') as fh:
                stats[asp] = json.load(fh)
        per = {asp: variant_inputs(a.variant, stats[asp], stats[asp]['tau'], {}, constrained) for asp in stats}
        act = set().union(*(per[x][0] for x in per))
        cost = {}
        for x in per:
            for g, c in per[x][1].items():
                cost[g] = cost.get(g, 0.0) + c
        truth = get_true_lineage(taxon, parents)
        rec = {'taxon': taxon, 'n_active': len(act), 'saved_status': saved['status'],
               'required': saved['required'], 'saved_category': saved['category']}
        if not act:
            rec.update(witness=[], chain=[], deepest=None)
        else:
            t0 = time.monotonic()
            adj, assign, _ = tc.solve_sat_ortools_taxon_assignment(
                {'genome': act}, only, never, hierarchy, cost, isa, dis, uni, time_limit=a.time_limit)
            rec['seconds'] = round(time.monotonic() - t0, 2)
            demoted = sorted(g for g, keep in adj['genome'].items() if not keep)
            rec['n_demoted'] = len(demoted)
            rec['demoted_cost'] = sum(cost.get(g, 0) for g in demoted)
            rec['reproduces_saved'] = (len(demoted) == saved.get('n_demoted') and
                                       abs(rec['demoted_cost'] - saved.get('demoted_cost', 0)) < 1e-6)
            rec['witness'] = sorted(assign)
            chain = sorted(t for t in assign if not t.startswith('NCBITaxon_Union_'))
            rec['chain'] = chain
            anc = {t: tc.get_all_ancestors_of_list(isa, [t]) for t in chain}
            deepest = [t for t in chain if set(chain) - {t} <= anc[t]]
            if chain and len(deepest) != 1:
                raise RuntimeError(f'Selected taxa are not one lineage chain: {taxon} {chain}')
            rec['deepest'] = deepest[0] if chain else None
        chain = rec['chain']
        rec['required_subset_of_witness'] = set(saved['required']) <= set(rec.get('witness', []))
        ranked = [t for t in chain if ranks.get(t) in RANKS_7]
        order = {r: i for i, r in enumerate(RANKS_7)}
        ranked.sort(key=lambda t: order[ranks[t]])
        rec['ranked_chain'] = [(t, ranks[t], names.get(t)) for t in ranked]
        rec['deepest_ranked'] = ranked[-1] if ranked else None
        rec['deepest_ranked_rank'] = ranks[ranked[-1]] if ranked else None
        req_ranked = sorted((t for t in saved['required'] if ranks.get(t) in RANKS_7), key=lambda t: order[ranks[t]])
        rec['required_deepest_ranked'] = req_ranked[-1] if req_ranked else None
        rec['required_deepest_rank'] = ranks[req_ranked[-1]] if req_ranked else None
        d = rec['deepest']
        rec['deepest_rank'] = ranks.get(d) if d else None
        rec['deepest_name'] = names.get(d) if d else None
        rec['deepest_in_true_lineage'] = (d in truth) if d else None
        rec['deepest_beyond_required'] = (d not in saved['required']) if d else None
        rec['deepest_ranked_beyond_required'] = (ranked[-1] != rec['required_deepest_ranked']) if ranked else None
        rec['chain_all_in_true_lineage'] = all(t in truth for t in chain)
        rec['fallback'] = not ranked
        done[taxon] = rec
        print(a.predictor, taxon, rec['deepest_name'], rec['deepest_rank'], 'inTrue', rec['deepest_in_true_lineage'],
              'beyondReq', rec['deepest_beyond_required'], 'repro', rec.get('reproduces_saved'), flush=True)
        out.write_text(json.dumps(done))
    out.write_text(json.dumps(done))


if __name__ == '__main__':
    main()
