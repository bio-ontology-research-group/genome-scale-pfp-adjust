# Reproducing the revised analysis

The October 2026 audit found that the earlier reproduction script copied tables and figures from the manuscript and silently omitted missing inputs. Those outputs are not accepted evidence for the revision.

The adjustment workflow consumes identified predictor inputs. Corrected local Seq-Sim and MLP-ESM2 generation is documented in [BASELINES.md](BASELINES.md); DeepGO-SE and SPROF-GO remain supplied pretrained outputs. The complete corrected numerical bundle and its public archive are still being prepared; the repository does not currently contain all raw prediction files.

## Required inputs

Use `revision.env.example` to identify an immutable code checkout, the prediction/annotation root, versioned ontology resources, the Python environment, a GAEF checkout, published FunTaxIS-lite downloads, complete source-proteome UniProt DAT files and the curated Complex Portal index. Do not replace an archived input with a current download and call it the same experiment.

The prediction root must contain `<predictor>_results/<aspect>/predictions/predictions_fold_*_taxon_<id>.tsv` and `annotations-go-basic/annots_taxon_<id>.tsv`. Exactly one prediction file is required for each task. The committed `revision_tasks.txt` specifies 1,764 tasks covering three GO aspects for the 50 entry-date organisms and 219 heldout organisms in their designated predictor comparisons.

The resource directory must contain `go-basic.obo`, `ncbitaxon.obo`, `go_taxon_constraints_extracted_obo.tsv`, `taxon_hierarchy.tsv` and `ncbitaxon_hierarchy.tsv`. Additional GAEF resource files are supplied in `data/constraints/` and checked by the input manifests. The audited GO snapshot is 2026-05-19; individual source-proteome files are identified by the reference-build manifest.

The official FunTaxIS-lite service downloads can be retrieved using `data_prep/fetch_funtaxis_constraints.py --help`. Preserve the original downloaded lists and manifest: a future service response need not reproduce the archived comparator.

## Rebuild reference labels before evaluation

The inherited annotation TSVs included propagation through relations such as
`regulates`. Safe repropagation cannot remove those extra labels. Rebuild the
references from the direct UniProt evidence, retaining the original corpus scope:

```bash
python analysis/rebuild_reference_annotations.py \
  --source-root /path/to/uniprot_reference_proteomes \
  --proteome-index data/uniprot_proteomes_ids.tsv \
  --original-annotations /path/to/original/annotations-go-basic \
  --go /path/to/versioned/go-basic.obo \
  --output /path/to/rebuilt-reference-annotations
```

Set `ANNOTATION_DIR` to that output directory. The driver requires its completed
`reference-manifest.json`; the aggregator verifies every reference checksum and
uses these labels for both evaluation and fixed-corpus information content.
The full reference corpus needs all source DATs recorded in that manifest, not
only the 50 DATs used in the separate evidence-sensitivity analysis.

Existing repair outputs can be reused because optimization does not consult
reference labels, but only when every evaluation protein has a stored score
vector in every condition. The aggregator checks that coverage and the original
and rebuilt reference hashes. An uncovered protein requires rerunning its task;
it must not silently receive a fabricated zero vector.

## Run the analysis

The canonical CPU analysis uses Python 3.10.20 and the packages pinned in
`requirements-revision.txt`. Install them in a dedicated environment; no GPU,
PyTorch, sequence aligner or upstream model training is needed for this workflow.
The broader `environment.yml` is for the older training stack and is not the
canonical revision environment. Recorded versions do not guarantee identical
wall-clock-bounded solver incumbents on different hardware.

```bash
python3.10 -m venv /your/project/revision-venv
/your/project/revision-venv/bin/pip install -r requirements-revision.txt
```

```bash
cp reproducibility/revision.env.example /your/project/revision.env
# Set the absolute input and output paths, environment and cluster account.
bash reproducibility/run_all_slurm.sh /your/project/revision.env
```

The driver submits a bounded worker array, evidence analysis, summary aggregation, lineage analysis, curated-membership validation and table generation. Dependencies use `afterok`; a failed task blocks the final generation job. It uses one CPU per process and respects the specified worker count. `AUDIT_MEMORY` controls worker memory (64G by default for the denser corrected baseline exports). `AUDIT_AGGREGATE_MEMORY` and `AUDIT_LINEAGE_MEMORY` control postprocessing memory (128G each by default). Dense MLP result metadata can exceed 1 GB compressed for a single organism; decoding the repeated score vectors needs substantially more memory than their on-disk size. These resource requests are recorded in the frozen run configuration. All scratch files go under the output directory. The source checkout and inputs must remain unchanged while jobs run; aggregation rejects mixed source manifests.

The workflow saves job IDs as it submits them. A scheduler rejection does not undo jobs already submitted: inspect `jobs.txt` and the queue before retrying. Use a new output directory for a new complete run. Failed or incomplete runs must not be reported as complete.

A single task can also be run with `analysis/revision_audit_run.py --help`. The strict entry-date protein mask applies to accuracy only; structural metrics retain all proteins. Known-taxon consistency and the GAEF pairwise compatibility screen are separate outputs. Stage 2 records feasible/optimal status and objective bounds, with no successful fallback on solver failure.

## Generate tables from a complete result bundle

```bash
python reproducibility/reproduce_all_tables_and_figures.py \
  --results /path/to/complete-results \
  --output /path/to/generated
```

Generation requires all 15 predictor/aspect summaries, `lineage.json`, `evidence-sensitivity.json` six `membership_<predictor>_<mode>.json` reports and `export-validation.json` covering all 300 serialized hard repairs. Result and prediction hashes must agree across those analyses. Missing files or inconsistent counts cause an error. TSV and LaTeX tables are produced from the same data, along with input/output checksums. The generator publishes its manifest last, after every artifact has been produced successfully. No manuscript PDF, table text or image is used as a numerical source.

| Artifact | Source |
|---|---|
| Main Table 1: counts and denominators | Native known-taxon Stage 1 summaries |
| Main Table 2: paired accuracy and known-taxon consistency | Native raw and Stage 1 summaries |
| S1: candidate/matched/eligible populations | Explicit scope counters |
| S2: raw, known and inferred-taxon accuracy | Native summaries |
| S3–S4: eight-threshold counts, consistency and accuracy | Complete threshold summaries, including MLP-ESM2 |
| S5: all seven lineage ranks | Conditional entailment checks and joint/aspect comparisons |
| S6: published FunTaxIS-lite and soft objectives | Comparator summaries |
| S7: evidence sensitivity | Complete-proteome annotation audit |
| S8: curated memberships | Exact accession mapping and paired Stage 1/2 outputs |
| S9: solver status and bounds | Per-organism CP-SAT records |
| S10: other structural measures before/after across thresholds | Complete threshold summaries |
| S11: GO identifier alignment | Per-input alias, merge and exclusion counts |
| Figure S1: illustrative cardinality repair | Executable TikZ source in `figures/repair-schematic.tex` |

The diagram is an explicitly labeled numerical example, not an empirical result. Compile it in the supplementary LaTeX document with TikZ. Historical metagenomic figures, stale numerical tables and the unsupported exact-reduction claim are not part of this revised evidence bundle.

## Regression checks

```bash
python -m pip install -r requirements-test.txt
mkdir -p "$HOME/tmp/genoadjust-tests"
export TMPDIR="$HOME/tmp/genoadjust-tests"
python data/constraints/verify_checksums.py
python tests/test_data_loaders.py
python tests/test_solver_smoke.py
python tests/test_golden_e2e.py
python tests/test_revision_contracts.py
python tests/test_evaluation_contracts.py
python tests/test_lineage_contracts.py
python tests/test_reproduction_workflow.py
python tests/test_table_generation.py
```

These tests cover the score contract, taxonomy unions, hierarchy and singleton constraints, evaluation membership, false-positive vocabulary and PR-curve endpoints. They do not establish biological correctness or substitute for the complete experiment and manuscript checks.

Stage 2 uses a 120-second search budget. A task that returns CP-SAT `UNKNOWN` without a feasible incumbent is retried once with 900 seconds for all its Stage 2 conditions. The accepted JSON records the actual budget; attempt logs are retained. Other errors, or a failed retry, stop the dependency chain.

The revision runner maps alternative GO IDs to current primary IDs before every comparison, merges duplicate scores by maximum, propagates scores to ancestors by maximum before every comparison, and records added/raised ancestor scores separately from repair flips. It also records dropped obsolete/unknown and out-of-aspect pairs. The lower-level TSV-only adjustment APIs expect IDs already aligned to the supplied hierarchy; use the revision runner when starting from predictor files using older GO vocabularies.
