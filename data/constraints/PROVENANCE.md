# Constraint-file provenance and versioning

These files are vendored snapshots of the background knowledge the two solver
stages enforce. They are pinned by SHA-256 in `checksums.sha256` and verified by
`verify_checksums.py` (also run in CI). Regenerating any file changes its hash;
update `checksums.sha256` and this table together when that happens.

| File | Source | GO / source release |
|------|--------|---------------------|
| `taxon_constraints.tsv` | GO `only_in_taxon` / `never_in_taxon` constraints, expanded over the GO hierarchy | GO 2025-10 |
| `go-taxon-groupings.owl` | GO taxon-grouping classes (taxon union/disjointness scaffolding) | GO 2025-10 |
| `ncbitaxon_with_disjointness.owl` | NCBITaxon with added `disjoint_from` / `union_of` axioms used by Stage 1 | NCBITaxon snapshot, 2025-10 |
| `heteromeric_complexes_2025_10.tsv` | Obligate heteromeric complex GO terms (GAEF set) | GO 2025_10 (stated in file header) |
| `protein_complexes.tsv` | Complex-membership table for Stage 2 complex coherence | GO 2025-10 |
| `essential_terms.tsv` | Essential-term list for the (future) completeness axis | GO 2025-10 |
| `has_part_relations.txt` | GO `has_part` relations for the (future) process-coherence axis | GO 2025-10 |
| `ec2go_v2025-03-16` | EC-number to GO mapping | GO/EC mapping 2025-03-16 (in filename) |
| `metacyc_GO_v2025-03-16_with_EC.tsv` | MetaCyc to GO mapping with EC numbers | 2025-03-16 (in filename) |

Notes:
- The two `*2025-03-16` mapping files carry the release in their filenames; the
  remaining constraint files were generated from the GO 2025-10 release used for
  the experiments in the paper.
- The Gene Ontology is released monthly; archived releases are at
  <http://release.geneontology.org/> (e.g. `2025-10-01`). NCBITaxon releases are
  at <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/>.
- To re-pin after an update: `sha256sum <files> > checksums.sha256` from this
  directory, then bump the release column above.
