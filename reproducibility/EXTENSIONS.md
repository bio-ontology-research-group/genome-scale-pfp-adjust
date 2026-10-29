# Revision extensions (2026-10-29)

These experiments extend the validated revision at commit `716bbbf`. They do not
replace its original benchmark results. Parameters below are fixed before the
new results are inspected; no test-set tuning or outcome-based case selection is
performed.

## Conservative complex repair

`complex_coherence.conservative.repair` accepts canonical, ancestor-closed scores,
a GO hierarchy, the operational complex vocabulary, and evidence maps keyed by
`(GO term, protein)` with nonempty source descriptions. These maps authorize
edits; the caller is responsible for the scientific validity of the evidence.
Positive evidence authorizes ancestor promotions. Negative evidence authorizes
descendant demotions. Demotions are forbidden unless explicitly authorized.
Protected input annotations and their ancestors cannot be removed. Conflicting
negative evidence and protection are rejected.

The existing margin-weighted objective gains a singleton slack penalty (default
1.0). Hierarchy remains hard. Every unresolved singleton is reported, including
whether no edit was authorized or optimization/coupled constraints prevented a
repair. A feasible time-limited incumbent is explicitly distinguished from a
proven optimum. All supplied eligible candidates are considered: this mode does
not impose the standard top-50 truncation. Optimality remains relative to the
specified evidence policy and integer objective.

The proteome experiments cover all 50 entry-date organisms for each of the three
predictors, using their validated known-taxon Stage 1 exports. The primary
`corroborated` policy requires both other predictors to assign the candidate
pair scores above their own native thresholds. This permits an unscored primary
candidate when there is corroborating computational support. No curated
membership enters this gate. The predictors are not statistically independent,
and agreement is not evidence that proteins physically assemble. No demotions
are allowed.

The `conservative` nonzero-score ablation authorizes a
promotion for an initially singleton complex only when that pair has a strictly
positive upstream score, and allows no demotions. Native thresholds and singleton
penalty 1.0 are fixed. This is **score support, not independent biological
validation**. Seq-Sim's serialization floor of 0.001 leaves few or no positive
inactive pairs at its native threshold. Report this limitation rather than
interpreting abstention as accurate partner identification. Exported scores use
the same reconstruction as the validated revision.

## Uncertainty across optimal GO repairs

`solve_sat_ortools_taxon_assignment(..., uncertainty=report, query_taxa=taxa)` first
solves the original Stage 1 objective. Only a proven optimum supports an
optimal-repair certificate. The model is then restricted to that exact integer
cost, its objective is cleared, and feasibility checks test the opposite of each
taxon value in the optimum witness. A taxon is required only when its negation is
proven infeasible. Feasible alternatives establish ambiguity; timeouts remain
undetermined. A taxon outside the model is not asserted. The existing return
values and default behavior remain unchanged.

`near_optimal_units` optionally expands the cost bound in integer objective units
(1,000,000 units per score-cost unit). This is an objective tolerance, not a
probability or calibrated confidence interval. The revision experiment uses zero
extra units and a 120-second budget per solve. It covers all five predictor/split
combinations and all CC, MF, BP and joint outcomes, comparing the selected-repair
conditional entailments to entailments across all optimal repairs. Accuracy and
resolution are reported separately; abstention must not count as a correct
resolved assignment.

## Targeted recovery experiment

Use the retained Complex Portal snapshot, not the merged organism–GO lookup:
retain individual records with exactly two distinct protein accessions, a term
in the operational complex vocabulary, and unique mappings of both proteins to
the full input proteome. Deduplicate identical (term, accession-pair) records.
Use each member in turn as the one supplied positive, keeping its partner hidden
from the repair input. Report exclusions and the exact qualifying case set.

This is a **controlled singleton-completion stress test**, not a natural missing
annotation benchmark. For each case, retain the supplied member at score
`max(0.9, tau+0.001)`. Compress every other protein's target score monotonically to
`tau * original_score`, so every other target is inactive under the strict
threshold. Propagate these scores to the target's ancestors. The original
Stage 1 score remains separate candidate evidence. Conservative repair allows
only candidates whose original score exceeds the native threshold, protects the
supplied member, and uses penalty 1.0. Crucially, the hidden member's identity is
never an argument to the function constructing repair inputs or running solvers.
It is used only afterward to score outputs.

Also run conservative repair using agreement of both other predictors as
candidate evidence, with the same protected supplied member. The hidden curated
membership is not used by this gate. Record whether the hidden member has this
corroboration separately from whether the primary predictor supports it.

Compare standard top-50 repair, demotion-only repair, both conservative policies,
and a highest-original-score eligible candidate rule for each policy
(identifier tie-breaking).
All methods use the complete source protein universe. Report hidden-member
recovery, removal of the supplied positive, unresolved cases, other proposals,
and the fraction of hidden members already supported by the original predictor.
Other proposals are **not validated negatives**. Zero protected-member removals
is a policy guarantee, not an estimated biological improvement. Results for
related terms, complexes and the two directions are correlated; do not treat
cases as independent biological replicates. Predictor training exposure remains
unresolved for supplied pretrained models.

## Running

`analysis/revision_extensions.py` has `lineage`, `conservative`, `corroborated`, and `recovery`
subcommands; use `--help` for paths and fixed parameter options. Each result
records source and resource hashes, input identities, solver status, and
case-level outcomes. Conservative runs also export prediction TSVs. Write new
results to a separate directory; never overwrite the validated revision.

For a new proteome outside the experiment layout, use the general CLI:

```bash
python pipeline/conservative_repair.py \
  --predictions stage1.tsv --go /path/to/go-basic.obo --threshold 0.3 \
  --promotion-evidence candidate-evidence.json \
  --protected-annotations trusted-input-annotations.json \
  --output conservative.tsv --audit conservative-audit.json
```

Evidence files are JSON arrays, for example
`[{"protein":"PROTEIN_ID","term":"GO:0032991","source":"traceable evidence identifier"}]`.
No evidence file means no edits of that kind are authorized. An optional
`--demotion-evidence` file authorizes removals; negative evidence is never inferred
from a missing prediction. The CLI canonicalizes identifiers and propagates CC
scores before repair, recording those preprocessing changes separately. Supply a
Stage 1 file when known-taxon consistency is required.

Recovery exclusion counters before deduplication count parsed GO-annotation rows;
mapping exclusions count deduplicated term/accession-pair candidates. Counters
are repeated for each predictor. Protein-pair cases can contain additional
non-protein constituents in the source assembly; the challenge concerns the two
indexed protein members. The retained qualifying source records were checked
for additional participant identifiers when interpreting this population.

The Stage 1 Python API accepts an output dictionary without changing its normal
return values:

```python
from taxon_consistency.adjust_ortools import solve_sat_ortools_taxon_assignment

uncertainty = {}
adjusted, chosen_taxa, flips = solve_sat_ortools_taxon_assignment(
    active_terms, only_in_taxon, never_in_taxon, go_parents, term_costs,
    taxon_parents, taxon_disjoint, taxon_unions,
    uncertainty=uncertainty, query_taxa=ranked_taxa, time_limit=120,
)
required_taxa = uncertainty["required"]
undetermined_taxa = uncertainty["undetermined"]
```

Here `active_terms` maps protein or aspect identifiers to active GO sets,
`term_costs` contains the Stage 1 genome-level demotion costs, and `ranked_taxa`
contains taxonomy identifiers using the resource's `NCBITaxon_` spelling.
`chosen_taxa` is a single compatible assignment; `required_taxa` contains the
certified conclusions. The report also retains possible alternatives, query
statuses, the exact integer optimum, and the objective scale.
