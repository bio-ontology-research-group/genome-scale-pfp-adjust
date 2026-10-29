"""Render the 2026-10-29 report; reject changed data before reusing its interpretation."""
import argparse
import hashlib
import json
from pathlib import Path

MODELS={'seq_sim_timeset':'Seq-Sim','deepgo-se':'DeepGO-SE','sprof-go':'SPROF-GO'}
POLICIES={'hard_standard':'Hard standard','hard_demotion':'Demotion-only','conservative':'Nonzero gate','corroborated':'Corroboration gate'}
RECOVERY={'standard':'Hard standard','demotion_only':'Demotion-only','conservative':'Primary-score conservative',
          'highest_score':'Primary-score greedy','conservative_corroborated':'Corroborated conservative',
          'highest_score_corroborated':'Corroborated greedy'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--summary',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args(); d=json.loads(args.summary.read_text())
    # Input paths can change when the review bundle is moved; scientific data cannot.
    scientific = {k: v for k, v in d.items() if k not in ('inputs', 'source_sha256')}
    digest = hashlib.sha256(json.dumps(scientific, sort_keys=True).encode()).hexdigest()
    if digest != 'ae974a73a134c263e5b93ee2f1cc8a36a27df8012c93cd8b4021394356e4ffd6':
        raise ValueError('Results changed: review the report interpretation before rendering')
    assert d['counts']['lineage_outcomes']==2352 and d['counts']['validated_exports']==300
    l=d['lineage']; q=sum(r['queries'] for r in l); c=sum(r['robust_claims'] for r in l)
    assert sum(r['cases_losing_claims'] for r in l)==0 and sum(r['undetermined'] for r in l)==0
    lines=['# Validated revision extension results','',
        'The conservative mode preserves supplied annotations and exposes unresolved cases. '
        'The optimal-repair analysis establishes that the reported required ranks are invariant across optimal GO repairs. '
        'The recovery experiment finds no improvement over the corresponding highest-score controls.','',
        '## Lineage conclusions across optimal repairs','',
        f'All 2,352 optimization problems reached proven optima. All {q:,} feasibility checks resolved, '
        f'and all {c:,} ranked assertions required by selected repairs remained required across every optimal repair. '
        'The reported lineage accuracy and resolution therefore remain unchanged. The correct joint ancestral-prefix '
        'counts are 41/50 for DeepGO-SE, 10/50 for entry-date Seq-Sim, 17/50 for SPROF-GO, '
        '25/219 for heldout Seq-Sim and 41/219 for MLP-ESM2. These are coarse lineage conclusions; '
        'the deepest correct joint rank is order.','',
        'This analysis removes dependence on the selected optimum for the reported required taxa. '
        'Incorrect conclusions remain consequences of the encoded inputs and objective. Certificates refer to '
        'the integer objective and represented taxonomy, with a 120-second budget per solve.','',
        '## Conservative proteome repair','',
        'All 300 conservative exports passed independent checks of input preservation, evidence restrictions, '
        'hierarchy, known taxonomy, score reconstruction and provenance. Corroboration requires both other '
        'predictors above their native thresholds; the nonzero ablation uses positive inactive primary scores. '
        'Both prohibit deletions in this experiment and penalize unresolved singletons with weight 1.','',
        '| Predictor | Policy | Added | Deleted | Unresolved / initial | Listed / assessable additions | Listed deletions |',
        '|---|---|---:|---:|---:|---:|---:|']
    for p in MODELS:
        for mode in POLICIES:
            r=next(r for r in d['complexes'] if r['predictor']==p and r['mode']==mode)
            lines.append(f"| {MODELS[p]} | {POLICIES[mode]} | {r['promotions']:,} | {r['demotions']:,} | {r['unresolved']:,} / {r['initial_singletons']:,} | {r['listed']} / {r['assessable']} | {r['listed_demotions']} |")
    cross=[r for r in d['complexes'] if r['mode']=='corroborated']
    unresolved=sum(r['unresolved'] for r in cross); initial=sum(r['initial_singletons'] for r in cross)
    added=sum(r['promotions'] for r in cross)
    lines+=['',f'Corroboration resolved {initial-unresolved:,} of {initial:,} initial singletons using {added:,} additions '
             f'and left {unresolved:,} unresolved. All 150 corroboration solves were optimal within their admitted '
             'candidate sets. For the nonzero ablation, 114/150 solves were optimal and 36 DeepGO-SE solves '
             'were feasible incumbents at the 120-second limit.','',
        'The zero-deletion result follows from the preservation policy. Seq-Sim and SPROF-GO have no admitted '
        'positive inactive candidates in the nonzero ablation. The DeepGO-SE nonzero gate allows many more '
        'additions than the hard repair. Corroboration is computational support, and agreement between '
        'predictors is not independent validation of physical assembly. A missing curated membership '
        'is not a validated negative.','',
        '## Controlled hidden-member recovery','',
        'The independently reconstructed experiment contains 292 directional cases per predictor, '
        '876 in total, from four organisms. One indexed protein member is supplied; its partner is withheld '
        'from the repair input. Other target scores are compressed to the inactive range while original '
        'predictions supply candidate evidence. Hidden membership is used only for scoring outcomes. '
        'Both directions and related complex terms are correlated cases.','',
        '| Policy | Seq-Sim recovered / 292 | DeepGO-SE recovered / 292 | SPROF-GO recovered / 292 | Supplied members deleted | Unresolved cases |',
        '|---|---:|---:|---:|---:|---:|']
    for mode,name in RECOVERY.items():
        rows=[next(r for r in d['recovery'] if r['predictor']==p and r['mode']==mode) for p in MODELS]
        values=' | '.join(f"{r['recovered']} / {r['cases']}" for r in rows)
        lines.append(f"| {name} | {values} | {sum(r['seed_removed'] for r in rows)} | {sum(r['unresolved'] for r in rows)} |")
    lines+=['',
        'The primary-score conservative policy recovered 122 hidden members, compared with 124 for its greedy '
        'control. Corroborated conservative repair and its greedy control each recovered 15. All these '
        'policies preserved supplied members. Standard hard repair recovered 136 and deleted 111 supplied '
        'members. The demotion-only result follows directly from requiring a hard repair of a constructed '
        'singleton while forbidding additions.','',
        'Primary scores supported 273 hidden members, and corroborating predictors supported 63. '
        'The additional evidence allowed SPROF-GO to recover 13 members where its primary-score gate '
        'recovered none; its corroborated greedy control recovered the same 13. The optimizer supplied '
        'no recovery advantage over these controls in the tested settings. Other proposed proteins are '
        'reported separately from recovery of the specified hidden member.','',
        'This is a controlled completion test with transformed scores. It does not establish prospective '
        'biological accuracy, and pretrained-model training exposure remains unresolved. The preservation '
        'policy and optimal-repair entailments are the defensible additions to the revision.','',
        '## Files and reproduction','',
        '- `summary.json`: complete summaries, input identities and validation links.',
        '- `provenance.json`: verified executed snapshots and post-execution validation changes.',
        '- `lineage-summary.tsv` and `lineage-ranks.tsv`: 20 groups and 16,464 ranked comparisons.',
        '- `complex-summary.tsv` and `complex-cases.tsv`: hard comparators and 300 conservative outputs.',
        '- `recovery-summary.tsv` and `recovery-cases.tsv`: six comparators and every recovery outcome.','',
        'The executed R1 lineage and R2 complex/recovery source snapshots, raw per-task results and independent '
        'validation reports are retained in the review bundle. The final repository also contains a general '
        'conservative-repair CLI, input-validation guards and report generators. The experiments use the '
        'previously validated Stage 1 inputs and the same curated-reference snapshot.','',
        'See [the extension protocol](../../../reproducibility/EXTENSIONS.md) for the evidence contracts, '
        'fixed parameters and execution commands. No gate or penalty was tuned against these recovery outcomes.','']
    args.output.write_text('\n'.join(lines))

if __name__=='__main__':
    main()
