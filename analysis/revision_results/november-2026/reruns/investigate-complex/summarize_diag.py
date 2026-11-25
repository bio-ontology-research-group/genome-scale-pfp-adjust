import gzip, json, sys, statistics as S
from collections import Counter
from pathlib import Path
root = Path(sys.argv[1]); out = {}
def med(x): return S.median(x) if x else None
for pdir in sorted(root.iterdir()):
    c = Counter(); lone_sc = []; dem_cost = []; prom_cost_alt = []; prom_rank = []; cp_rank = []; uni_rank = []; ratio = []
    ex_dem = []
    for f in sorted(pdir.glob('*.json.gz')):
        d = json.load(gzip.open(f, 'rt'))
        c['taxa'] += 1
        for r in d['records']:
            c['singletons'] += 1
            c['no_pos_inactive'] += r['n_pos_inactive'] == 0
            c['top50_all_zero'] += r['n_top50_pos'] == 0
            c['tie_at_min_cost_gt1'] += r['n_cand_at_min_cost'] > 1
            if not r['std_lone_kept']:
                c['std_demoted'] += 1
                dem_cost.append(r['dem_cost']); prom_cost_alt.append(r['best_prom_cost'])
                lone_sc.append(r['lone_score'])
                if r['best_prom_cost'] is not None:
                    c['std_dem_cheaper_than_best_prom'] += r['dem_cost'] <= r['best_prom_cost'] + 1e-9
                    ratio.append(r['dem_cost'] / r['best_prom_cost'])
                c['std_dem_lone_uni'] += r['uni']['lone_class'] is not None
                c['std_dem_lone_uni_curated'] += r['uni']['lone_class'] == 'curated'
                if 'cp' in r:
                    c['std_dem_cp_key'] += 1; c['std_dem_cp_lone_listed'] += r['cp']['lone_listed']
                    sub = [s for s in r['cp']['subunits'] if not s['is_lone']]
                    c['std_dem_cp_partner_in_proteome'] += bool(sub)
                    c['std_dem_cp_partner_in_top50'] += any(s['in_top50'] for s in sub)
                    for s in sub:
                        if len(ex_dem) < 15: ex_dem.append([f.stem, r['term'], r['lone'], round(r['lone_score'],3), round(r['dem_cost'],3), s['entry'], round(s['score'],3), s['rank'], round(s['prom_cost'],3), round(r['best_prom_cost'],3)])
            if r['std_promoted']:
                c['std_promoted_X'] += 1
                prom_rank += r['std_promoted_rank']
                c['std_prom_pairs'] += len(r['std_promoted'])
                c['std_prom_score0'] += sum(s == 0 for s in r['std_promoted_score'])
                c['std_prom_at_min_cost'] += sum(abs(pc - r['best_prom_cost']) < 1e-9 for pc in r['std_promoted_cost'])
            if 'cp' in r:
                c['cp_key'] += 1
                c['cp_lone_listed'] += r['cp']['lone_listed']
                c['cp_listed_total'] += r['cp']['n_listed']; c['cp_listed_mapped'] += r['cp']['n_mapped']
                sub = [s for s in r['cp']['subunits'] if not s['is_lone']]
                c['cp_has_partner_in_proteome'] += bool(sub)
                c['cp_partner_any_in_top50'] += any(s['in_top50'] for s in sub)
                c['cp_partner_any_score_pos'] += any(s['score'] > 0 for s in sub)
                c['cp_partners'] += len(sub); c['cp_partners_score0'] += sum(s['score'] == 0 for s in sub)
                c['cp_partners_in_top50'] += sum(s['in_top50'] for s in sub)
                cp_rank += [s['rank'] for s in sub]
                if sub and r['best_prom_cost'] is not None:
                    c['cp_partner_cheapest_le_best'] += min(s['prom_cost'] for s in sub) <= r['best_prom_cost'] + 1e-9
                c['cp_std_promoted_listed'] += bool(set(r['std_promoted']) & {s['entry'] for s in sub})
                c['cp_std_promoted_any'] += bool(r['std_promoted'])
            um = r['uni']['members']
            c['uni_has_partner'] += bool(um)
            c['uni_has_partner_curated'] += any(m['cls'] == 'curated' for m in um)
            c['uni_partner_in_top50'] += any(m['in_top50'] for m in um)
            c['uni_partner_score_pos'] += any(m['score'] > 0 for m in um)
            uni_rank += [m['rank'] for m in um]
            c['uni_lone_member'] += r['uni']['lone_class'] is not None
            c['uni_std_promoted_member'] += bool(set(r['std_promoted']) & {m['entry'] for m in um})
    out[pdir.name] = dict(counts=dict(c), median_lone_score_demoted=med(lone_sc), median_dem_cost=med(dem_cost),
        median_best_prom_cost_when_demoted=med([x for x in prom_cost_alt if x is not None]), median_dem_over_prom=med(ratio),
        median_rank_promoted=med(prom_rank), median_rank_cp_partner=med(cp_rank), median_rank_uni_partner=med(uni_rank),
        frac_uni_partner_rank_le50=sum(x <= 50 for x in uni_rank)/len(uni_rank) if uni_rank else None,
        cp_demotion_examples=ex_dem)
print(json.dumps(out, indent=1))
