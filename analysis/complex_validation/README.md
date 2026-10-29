# Curated complex-membership checks

Compare a Stage 2 output with its corresponding Stage 1 input. The analysis counts all edited protein--term pairs, but assesses curated membership only for terms in the same operational complex vocabulary used by the repair. Complex Portal also supplies process, function and location annotations; these must not be counted as complex-membership references.

UniProt entry names require an explicit accession map. A changed pair is assessable only when its organism--complex reference and protein accession are available. Report the total changes, reference coverage, listed-subunit promotions, listed-subunit demotions and unresolved accession counts. An unlisted protein is not a validated non-member, and a listed member does not establish that two promoted proteins physically interact. Zero coverage means precision is unmeasured.

The accepted revision uses the retained Complex Portal download snapshot. CORUM is supported by the index builder but is not included in that snapshot. No QuickGO-derived labels are mixed into this independent comparison.

```bash
python analysis/complex_validation/build_subunit_index.py \
  --complexportal_dir /path/to/complexportal \
  --output_file /path/to/subunit_index.json

python analysis/complex_validation/validate_stage2.py \
  --stage1_dir /path/to/results/deepgo-se/cc/stage1 \
  --stage2_dir /path/to/results/deepgo-se/cc/stage2_standard \
  --subunit_index /path/to/subunit_index.json \
  --protein-id-map /path/to/results/protein-id-map.json \
  --go /path/to/exact/go-basic.obo \
  --threshold 0.3 --label 'DeepGO-SE standard' \
  --output_file /path/to/membership.json
```

The other two entry-date predictors use threshold 0.001. Run both standard and demotion-only comparisons for all three predictors. The result JSON records the ontology, exclusion-list, index, accession-map, script and paired prediction hashes. The table generator checks these paired files and edit totals against independent export validation. `reproducibility/run_all_slurm.sh` schedules the complete analysis with the required dependencies.
