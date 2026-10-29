# Audited revision analysis

The October 2026 audit invalidated the earlier numerical tables and several analysis claims. The historical files under `revision_results/` are retained for traceability; they are not the accepted results of the corrected workflow. In particular, the earlier exact sparse-reduction claim, locally approximated FunTaxIS-lite comparison, heuristic SBR label and unmasked entry-date accuracy must not be reused.

The current workflow and table-to-source mapping are in [reproducibility/README.md](../reproducibility/README.md). It consumes identified predictor inputs. Local baseline regeneration is documented in [BASELINES.md](../reproducibility/BASELINES.md); supplied pretrained predictions are distinguished from regenerated local baselines.

| Analysis | Current implementation |
|---|---|
| Canonical input alignment, threshold sweep, hard/soft repair, published filtering | `revision_audit_run.py` |
| Explicit test populations and fixed-corpus accuracy | `evaluation_scope.py`, `aggregate_revision_audit.py` |
| Seven-rank lineage entailment and cross-aspect agreement | `audit_lineage_results.py` |
| Complete-proteome evidence sensitivity and accession mapping | `audit_annotation_evidence.py` |
| Paired curated complex-membership coverage | `complex_validation/validate_stage2.py` |
| Independent serialized hard-repair checks | `validate_repair_exports.py` |
| Numerical tables and illustrative figure source | `../reproducibility/reproduce_all_tables_and_figures.py` |

A complete run includes 1,764 model/aspect/taxon tasks, all 15 aggregate groups, lineage and evidence analyses, six curated-membership comparisons and 300 independently checked hard-repair files. Missing or inconsistent results block final generation. Numerical acceptance and archival verification remain pending until that complete run passes.
