# Empty Quarter IBEX genome demo

*2026-07-01T05:37:34Z by Showboat 0.6.1*
<!-- showboat-id: 3cc37504-24d0-43e9-9e08-2196ae8dadaf -->

This demo runs GenoAdjust on a real Empty Quarter genome annotation from IBEX: rh04 / 63_rh04, identified by ANI as Bacillus spizizenii (NCBITaxon_96241). It fetches the PGAP GFF and ANI report from IBEX, converts PGAP's GO-bearing CDS annotations into the repository prediction format, runs taxon consistency, then runs optimized complex coherence with top_k=2. Raw genome sequence is not stored in the repository; generated inputs, outputs, and solver logs are written under the ignored output/ directory.

Requirements: an ssh/scp host alias named ibex with access to /ibex/scratch/projects/c2014/EmptyQuarter_Data, plus uvx so the pinned OR-Tools dependency can be provided without modifying the system Python.

```bash
uvx --with ortools==9.15.6755 --with tqdm python demos/empty_quarter_ibex/run_demo.py
```

```output
Source genome
  sample: rh04 / 63_rh04
  IBEX GFF: /ibex/scratch/projects/c2014/EmptyQuarter_Data/cultures/assemblies/assemblies_tiannyu_2025/pgap_annotations/63_rh04/annot.gff
  ANI best match: Bacillus spizizenii (taxid = 96241, rank = species, lineage = Bacteria; Bacillati; Bacillota; Bacilli; Bacillales; Bacillaceae; Bacillus)
  ANI status: CONFIRMED (HIGH)

Converted PGAP GO predictions
  CDS records: 4193
  CDS records with GO terms: 1999
  proteins with GO predictions: 1989
  distinct GO terms: 1280
  above-threshold annotations before adjustment: 4943

Stage 1: taxon consistency
  above-threshold annotations after Stage 1: 4932
  demotions: 11
  promotions: 0
  GO:0044423 (virion component): 7
  GO:0055051 (ATP-binding cassette (ABC) transporter complex, integrated substrate binding): 2
  GO:0000311 (plastid large ribosomal subunit): 1
  GO:0045087 (innate immune response): 1

Stage 2: complex coherence
  top_k: 2
  above-threshold annotations after Stage 2: 4911
  demotions: 21
  promotions: 0
  GO:0000015 (phosphopyruvate hydratase complex): 1
  GO:0005854 (nascent polypeptide-associated complex): 1
  GO:0005948 (acetolactate synthase complex): 1
  GO:0005951 (carbamoyl-phosphate synthase complex): 1
  GO:0005960 (glycine cleavage complex): 1
  GO:0005971 (ribonucleoside-diphosphate reductase complex): 1
  GO:0009318 (exodeoxyribonuclease VII complex): 1
  GO:0009320 (phosphoribosylaminoimidazole carboxylase complex): 1

Generated files
  predictions: demos/empty_quarter_ibex/output/predictions/predictions_fold_00_taxon_96241.tsv
  taxon-adjusted: demos/empty_quarter_ibex/output/optimized/taxon_fold_00_taxon_96241.tsv
  final adjusted: demos/empty_quarter_ibex/output/optimized/optimized_fold_00_taxon_96241.tsv
  raw solver logs: demos/empty_quarter_ibex/output/logs
```

The fixed score of 0.9 is only a representation layer: PGAP's GFF records contain discrete GO annotations rather than probabilities, while GenoAdjust expects GO|score entries. The demo therefore tests the repair code on real genome-derived GO predictions without committing the underlying IBEX data.
