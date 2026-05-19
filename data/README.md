# Bundled data

These small constraint files are committed to the repo; the larger ones are
fetched by `../data_prep/download_data.sh`.

## Bundled (committed)

| File | Source |
|---|---|
| `go-computed-taxon-constraints.obo` | Gene Ontology (release 2025-10) |
| `go-taxon-groupings.obo` | Gene Ontology |
| `go_hierarchy.tsv` | Parsed `is_a` edges from `go-basic.obo` |
| `go_hierarchy_{bp,cc,mf}.tsv` | Per-ontology splits |
| `go_taxon_constraints_updated.tsv` | Flattened `only_in_taxon`/`never_in_taxon` constraints |
| `go_taxon_constraints_extracted_obo.tsv` | Same, alternate extraction (kept for `analysis/`) |
| `taxon_hierarchy.tsv` | Curated NCBI Taxonomy subset for taxon unions |
| `uniprot_proteomes_ids.tsv` | TaxonID <-> ProteomeID mapping (regenerated via `data_prep/obtain_proteomes_ids.py`) |
| `../complex_coherence/protein_complexes.tsv` | 1,982 obligate heteromeric complexes from GAEF |
| `constraints/essential_terms.tsv` | Vendored from upstream GAEF |
| `constraints/has_part_relations.txt` | Vendored from upstream GAEF |
| `constraints/ec2go_v2025-03-16` | Vendored from upstream GAEF |
| `constraints/metacyc_GO_v2025-03-16_with_EC.tsv` | Vendored from upstream GAEF |
| `constraints/ncbitaxon_with_disjointness.owl` | Vendored from upstream GAEF |
| `constraints/go-taxon-groupings.owl` | Vendored from upstream GAEF |
| `constraints/protein_complexes.tsv` | Vendored from upstream GAEF (duplicates `../complex_coherence/protein_complexes.tsv`) |
| `constraints/taxon_constraints.tsv` | Regenerated against GO 2025-10 via GAEF's `generate_constraints/extract_constraints.groovy` (diverges from upstream's committed copy) |

## Downloaded (not committed)

| File | Source | How to regenerate |
|---|---|---|
| `go-basic.obo` (~31 MB) | `release.geneontology.org` | `../data_prep/download_data.sh` |
| `ncbitaxon.obo` (~700 MB) | `purl.obolibrary.org/obo/ncbitaxon.obo` | `../data_prep/download_data.sh` |
| `ncbitaxon_hierarchy.tsv` (~93 MB) | Parsed from `ncbitaxon.obo` | `../data_prep/download_data.sh` (calls `data_prep/parse_ncbitaxon.py`) |

## Regenerating from primary sources

```
./data_prep/download_data.sh                                                    # fetches go-basic.obo + ncbitaxon.obo (also runs parse_ncbitaxon)
python data_prep/parse_go_hierarchy.py    --go-obo-file data/go-basic.obo --output-file data/go_hierarchy.tsv
python data_prep/parse_taxon_hierarchy.py --go-taxon-groupings-file data/go-taxon-groupings.obo --output-file data/taxon_hierarchy.tsv
python data_prep/parse_ncbitaxon.py       --ncbitaxon-obo-file data/ncbitaxon.obo --output-file data/ncbitaxon_hierarchy.tsv
```
