# Stage 2 biological-validity validation against curated complex subunit lists

## What this checks

Stage 2 (complex coherence) repairs each "singleton complex" violation by
either **demoting** the term from the lone protein or **promoting** it onto
an additional protein. The choice is made by the solver on cost (score
margin), not on biology. So the headline "100% complex coherence after Stage
2" is a *structural* property — it does not say whether the promoted
partners are actually subunits of the complex.

This analysis closes that loop by comparing every Stage-2 promotion and
demotion against externally curated subunit lists:

- **Promotion precision**:  fraction of `(organism, complex GO term, protein)`
  promotions where the protein is a curated subunit of the complex.
- **Demotion characterisation**: fraction of demoted singletons that are
  themselves curated subunits (suggesting we removed a real annotation
  because its partners are simply not yet curated in the organism) versus
  non-subunits (likely upstream false positives that we correctly killed).

## Data sources

| Source | Coverage | Why we use it |
|---|---|---|
| EBI **Complex Portal** | E. coli, yeast, *D. melanogaster*, *C. elegans*, human, mouse, rat, *A. thaliana*, more | Multi-organism, curated, free download. Primary source. |
| **CORUM** core complexes | Mammalian (mostly human/mouse/rat) | Gold-standard mammalian set; secondary source. |
| (Optional) QuickGO experimental annotations | Any organism with curated GO annotations | Coverage extension for organisms that neither ComplexPortal nor CORUM cover. Uses the same evidence-code filter (EXP/IDA/IPI/IMP/IGI) as the paper's ground truth. |

The 50 timeset organisms include several that ComplexPortal covers
(E. coli K-12 strains, *S. cerevisiae*, *D. melanogaster*, …). Coverage
will be partial — that is the point: report `precision_on_covered` and
state the coverage fraction.

## Coverage caveat

A promotion landing on a non-curated protein is not necessarily wrong; the
subunit may not yet have been characterised. Treat the reported precision
as a **lower bound** on biological accuracy, and discuss this explicitly in
the paper.

## Workflow

```
1.  download_complex_dbs.sh   # ComplexPortal TSVs + CORUM coreComplexes
2.  build_subunit_index.py    # -> subunit_index.json  ({"<taxon>::<GO term>": [UniProts]})
3.  validate_stage2.py        # diff predictions vs optimized, classify against index
```

Outputs into `${DATA_DIR}/complex_validation/`:

- `subunit_index.json` — unified `(taxon, GO complex term) -> [UniProt subunits]` lookup
- `validated_pairs.tsv` — one row per Stage-2 action: `taxon, GO term, protein, action, curated_status`
- `summary.tsv` — aggregate precision/recall of promotions and demotions
- `per_organism.tsv`, `per_complex.tsv` — breakdowns

## How to invoke

```bash
source ../../config.env

bash download_complex_dbs.sh

python build_subunit_index.py \
    --complexportal_dir ${DATA_DIR}/complex_validation/complexportal \
    --corum_file       ${DATA_DIR}/complex_validation/corum/coreComplexes.txt \
    --output_file      ${DATA_DIR}/complex_validation/subunit_index.json

# Example: validate the Seq-Sim timeset Stage-2 outputs (CC only).
python validate_stage2.py \
    --predictions_dir  ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions \
    --optimized_dir    ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/optimized \
    --subunit_index    ${DATA_DIR}/complex_validation/subunit_index.json \
    --complexes_file   ../../complex_coherence/protein_complexes.tsv \
    --threshold        0.3 \
    --output_dir       ${DATA_DIR}/complex_validation/seq_sim_timeset_cc
```

Run once per predictor that has a Stage-2 output directory (Seq-Sim,
DeepGO-SE, SPROF-GO on CC for the timeset split).

## What goes into the paper

A single new paragraph in §3.6 reporting:

- For each predictor: `N_promoted_covered / N_promoted_total`,
  promotion precision, demotion characterisation, with a 95% binomial CI.
- A short Limitations sentence acknowledging the
  curation-completeness floor.
