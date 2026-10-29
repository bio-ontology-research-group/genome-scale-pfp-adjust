# Dataset splits

The revision preserves the committed organism splits and recorded training/validation protein membership. Corrected local baseline labels and vocabularies are rebuilt from source evidence as documented in `reproducibility/BASELINES.md`; pretrained DeepGO-SE and SPROF-GO prediction files remain fixed. The reference index contains 569 taxa; available annotation files cover 568 taxa.

## Entry-date split (`timeset/`)

`proteins_by_date_23-MAY-2024_filtered.tsv` lists 169 candidate Swiss-Prot entries integrated after 2024-05-23 in the 50 taxa in `test_organisms.txt`. The date is entry integration, not first sequence availability or first experimental annotation. This criterion does not establish exclusion from a pretrained model's training corpus.

Accuracy evaluation intersects exact entry names with the reference annotations and requires a non-root label in the evaluated GO aspect. Structural evaluation retains all proteins in each of the 50 prediction files. The generated scope table reports matched and eligible populations separately. For information content, candidate entries are excluded from the reference corpus.

The original extraction utility is `data_prep/extract_proteins_by_date.py`; rerunning it against a newer UniProt release creates a different input set and does not reproduce the archived experiment.

## Heldout split (`heldout/`)

`test_organisms.txt` contains 219 test taxa and `train_organisms.txt` contains 350 training taxa. The annotation file for training taxon 6645 is unavailable, so the fixed information-content corpus uses the other 349 training taxa. Every designated test task has its required prediction and annotation file.

`split_info.json` describes the original split construction. The revision does not independently establish a universal pairwise sequence-identity exclusion or freedom from pretrained-model leakage. Preserve the exact split files, supplied predictions and annotation snapshot when reproducing the reported results.
