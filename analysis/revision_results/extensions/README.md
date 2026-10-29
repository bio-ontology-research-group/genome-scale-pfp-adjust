# Validated revision extension results

The conservative mode preserves supplied annotations and exposes unresolved cases. The optimal-repair analysis establishes that the reported required ranks are invariant across optimal GO repairs. The recovery experiment finds no improvement over the corresponding highest-score controls.

## Lineage conclusions across optimal repairs

All 2,352 optimization problems reached proven optima. All 96,556 feasibility checks resolved, and all 4,890 ranked assertions required by selected repairs remained required across every optimal repair. The reported lineage accuracy and resolution therefore remain unchanged. The correct joint ancestral-prefix counts are 41/50 for DeepGO-SE, 10/50 for entry-date Seq-Sim, 17/50 for SPROF-GO, 25/219 for heldout Seq-Sim and 41/219 for MLP-ESM2. These are coarse lineage conclusions; the deepest correct joint rank is order.

This analysis removes dependence on the selected optimum for the reported required taxa. Incorrect conclusions remain consequences of the encoded inputs and objective. Certificates refer to the integer objective and represented taxonomy, with a 120-second budget per solve.

## Conservative proteome repair

All 300 conservative exports passed independent checks of input preservation, evidence restrictions, hierarchy, known taxonomy, score reconstruction and provenance. Corroboration requires both other predictors above their native thresholds; the nonzero ablation uses positive inactive primary scores. Both prohibit deletions in this experiment and penalize unresolved singletons with weight 1.

| Predictor | Policy | Added | Deleted | Unresolved / initial | Listed / assessable additions | Listed deletions |
|---|---|---:|---:|---:|---:|---:|
| Seq-Sim | Hard standard | 5,919 | 0 | 0 / 2,904 | 0 / 30 | 0 |
| Seq-Sim | Demotion-only | 0 | 2,949 | 0 / 2,904 | 0 / 0 | 8 |
| Seq-Sim | Nonzero gate | 0 | 0 | 2,904 / 2,904 | 0 / 0 | 0 |
| Seq-Sim | Corroboration gate | 13 | 0 | 2,901 / 2,904 | 0 / 0 | 0 |
| DeepGO-SE | Hard standard | 1,024 | 1,405 | 0 / 2,017 | 7 / 16 | 21 |
| DeepGO-SE | Demotion-only | 0 | 2,035 | 0 / 2,017 | 0 / 0 | 30 |
| DeepGO-SE | Nonzero gate | 4,242 | 0 | 22 / 2,017 | 17 / 58 | 0 |
| DeepGO-SE | Corroboration gate | 790 | 0 | 1,580 / 2,017 | 2 / 7 | 0 |
| SPROF-GO | Hard standard | 262 | 0 | 0 / 143 | 0 / 0 | 0 |
| SPROF-GO | Demotion-only | 0 | 143 | 0 / 143 | 0 / 0 | 0 |
| SPROF-GO | Nonzero gate | 0 | 0 | 143 / 143 | 0 / 0 | 0 |
| SPROF-GO | Corroboration gate | 0 | 0 | 143 / 143 | 0 / 0 | 0 |

Corroboration resolved 440 of 5,064 initial singletons using 803 additions and left 4,624 unresolved. All 150 corroboration solves were optimal within their admitted candidate sets. For the nonzero ablation, 114/150 solves were optimal and 36 DeepGO-SE solves were feasible incumbents at the 120-second limit.

The zero-deletion result follows from the preservation policy. Seq-Sim and SPROF-GO have no admitted positive inactive candidates in the nonzero ablation. The DeepGO-SE nonzero gate allows many more additions than the hard repair. Corroboration is computational support, and agreement between predictors is not independent validation of physical assembly. A missing curated membership is not a validated negative.

## Controlled hidden-member recovery

The independently reconstructed experiment contains 292 directional cases per predictor, 876 in total, from four organisms. One indexed protein member is supplied; its partner is withheld from the repair input. Other target scores are compressed to the inactive range while original predictions supply candidate evidence. Hidden membership is used only for scoring outcomes. Both directions and related complex terms are correlated cases.

| Policy | Seq-Sim recovered / 292 | DeepGO-SE recovered / 292 | SPROF-GO recovered / 292 | Supplied members deleted | Unresolved cases |
|---|---:|---:|---:|---:|---:|
| Hard standard | 104 / 292 | 32 / 292 | 0 / 292 | 111 | 0 |
| Demotion-only | 0 / 292 | 0 / 292 | 0 / 292 | 876 | 0 |
| Primary-score conservative | 104 / 292 | 18 / 292 | 0 / 292 | 0 | 485 |
| Primary-score greedy | 104 / 292 | 20 / 292 | 0 / 292 | 0 | 479 |
| Corroborated conservative | 1 / 292 | 1 / 292 | 13 / 292 | 0 | 671 |
| Corroborated greedy | 1 / 292 | 1 / 292 | 13 / 292 | 0 | 669 |

The primary-score conservative policy recovered 122 hidden members, compared with 124 for its greedy control. Corroborated conservative repair and its greedy control each recovered 15. All these policies preserved supplied members. Standard hard repair recovered 136 and deleted 111 supplied members. The demotion-only result follows directly from requiring a hard repair of a constructed singleton while forbidding additions.

Primary scores supported 273 hidden members, and corroborating predictors supported 63. The additional evidence allowed SPROF-GO to recover 13 members where its primary-score gate recovered none; its corroborated greedy control recovered the same 13. The optimizer supplied no recovery advantage over these controls in the tested settings. Other proposed proteins are reported separately from recovery of the specified hidden member.

This is a controlled completion test with transformed scores. It does not establish prospective biological accuracy, and pretrained-model training exposure remains unresolved. The preservation policy and optimal-repair entailments are the defensible additions to the revision.

## Files and reproduction

- `summary.json`: complete summaries, input identities and validation links.
- `provenance.json`: verified executed snapshots and post-execution validation changes.
- `lineage-summary.tsv` and `lineage-ranks.tsv`: 20 groups and 16,464 ranked comparisons.
- `complex-summary.tsv` and `complex-cases.tsv`: hard comparators and 300 conservative outputs.
- `recovery-summary.tsv` and `recovery-cases.tsv`: six comparators and every recovery outcome.

The executed R1 lineage and R2 complex/recovery source snapshots, raw per-task results and independent validation reports are retained in the review bundle. The final repository also contains a general conservative-repair CLI, input-validation guards and report generators. The experiments use the previously validated Stage 1 inputs and the same curated-reference snapshot.

See [the extension protocol](../../../reproducibility/EXTENSIONS.md) for the evidence contracts, fixed parameters and execution commands. No gate or penalty was tuned against these recovery outcomes.
