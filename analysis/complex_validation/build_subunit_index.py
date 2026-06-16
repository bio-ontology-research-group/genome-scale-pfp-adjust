"""
Build a unified curated-subunit lookup used by validate_stage2.py.

For each row in ComplexPortal / CORUM, we extract:
  - taxon ID
  - the GO complex term(s) the row is annotated with
  - the UniProt accessions of the subunits

and produce a JSON map keyed by "<taxon_id>::<GO term>" whose value is the
sorted list of curated UniProt subunits.

Coverage caveats are intentional: a missing entry means "no curated
reference for this (organism, complex)", not "no subunits exist". This is
why validate_stage2.py reports precision on the *covered* subset and reports
the coverage fraction alongside.
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict

# UniProt accession pattern (https://www.uniprot.org/help/accession_numbers).
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
GO_RE = re.compile(r"GO:\d{7}")


# ---------------------------------------------------------------------------
# ComplexPortal parsing
# ---------------------------------------------------------------------------
def _split_pipe_with_parens(field):
    """Split a ComplexPortal pipe-separated field, preserving stoichiometry parens."""
    if not field or field == "-":
        return []
    return [tok.strip() for tok in field.split("|") if tok.strip()]


def _extract_uniprot_from_token(tok):
    """ComplexPortal subunit tokens look like 'P12345(1)' or 'CHEBI:1234(1)' or 'EBI-12345(1)'.
    Keep only UniProt-shaped accessions."""
    acc = tok.split("(")[0].strip()
    return acc if UNIPROT_RE.match(acc) else None


def parse_complexportal_file(path):
    """Yield (taxon, go_term, frozenset(subunits)) entries from one complextab TSV."""
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        # ComplexPortal column names are stable across releases; if the layout
        # changes, fix the mapping here once.
        try:
            taxonomy_col = header.index("Taxonomy identifier")
            subunit_col = header.index(
                "Identifiers (and stoichiometry) of molecules in complex"
            )
            go_col = header.index("Go Annotations")
        except ValueError:
            sys.stderr.write(
                f"  [warn] {os.path.basename(path)}: unexpected header, skipping\n"
            )
            return

        for line in f:
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(taxonomy_col, subunit_col, go_col):
                continue
            taxon = row[taxonomy_col].split("(")[0].strip()
            if not taxon:
                continue
            subunits = set()
            for tok in _split_pipe_with_parens(row[subunit_col]):
                acc = _extract_uniprot_from_token(tok)
                if acc:
                    subunits.add(acc)
            if not subunits:
                continue
            go_terms = set(GO_RE.findall(row[go_col]))
            for go in go_terms:
                yield taxon, go, frozenset(subunits)


# ---------------------------------------------------------------------------
# CORUM parsing
# ---------------------------------------------------------------------------
CORUM_ORGANISM_TO_TAXON = {
    "Human": "9606",
    "Mouse": "10090",
    "Rat": "10116",
    "Pig": "9823",
    "Bovine": "9913",
    "Hamster": "10029",
    "Dog": "9615",
    "Rabbit": "9986",
    "Chicken": "9031",
}


def parse_corum_file(path):
    """Yield (taxon, go_term, frozenset(subunits)) entries from CORUM coreComplexes.txt."""
    # CORUM ships latin-1 in older releases; utf-8 in newer. Try both.
    try:
        f = open(path, encoding="utf-8")
        header = f.readline().rstrip("\n").split("\t")
    except UnicodeDecodeError:
        f = open(path, encoding="latin-1")
        header = f.readline().rstrip("\n").split("\t")

    # Column names have varied across CORUM releases. Match by best-known names.
    def find(col_candidates):
        for c in col_candidates:
            if c in header:
                return header.index(c)
        return None

    organism_col = find(["Organism", "organism"])
    subunit_col = find(
        [
            "subunits(UniProt IDs)",
            "subunits (UniProt IDs)",
            "Subunits (UniProt IDs)",
            "subunit(UniProtKB)",
        ]
    )
    go_col = find(["GO ID", "GO id", "FunCat ID"])  # FunCat is a fallback only if GO ID is missing
    if None in (organism_col, subunit_col, go_col):
        sys.stderr.write(
            f"  [warn] CORUM header layout unrecognised: {header}\n"
            "        Update CORUM column mapping in build_subunit_index.py.\n"
        )
        f.close()
        return

    for line in f:
        row = line.rstrip("\n").split("\t")
        if len(row) <= max(organism_col, subunit_col, go_col):
            continue
        taxon = CORUM_ORGANISM_TO_TAXON.get(row[organism_col].strip())
        if not taxon:
            continue
        subunits = set()
        # CORUM uses ';' or ',' as subunit separator.
        for s in re.split(r"[;,]", row[subunit_col]):
            s = s.strip()
            if UNIPROT_RE.match(s):
                subunits.add(s)
        if not subunits:
            continue
        for go in GO_RE.findall(row[go_col]):
            yield taxon, go, frozenset(subunits)
    f.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--complexportal_dir",
        required=False,
        help="Directory of ComplexPortal *.tsv files",
    )
    ap.add_argument(
        "--corum_file",
        required=False,
        help="Path to CORUM coreComplexes.txt",
    )
    ap.add_argument("--output_file", required=True, help="Output JSON file")
    args = ap.parse_args()

    index = defaultdict(set)
    sources = defaultdict(lambda: defaultdict(int))  # source -> {"entries": N, "taxa": set}

    if args.complexportal_dir and os.path.isdir(args.complexportal_dir):
        for fn in sorted(os.listdir(args.complexportal_dir)):
            if not fn.endswith(".tsv"):
                continue
            path = os.path.join(args.complexportal_dir, fn)
            for taxon, go, subs in parse_complexportal_file(path):
                index[f"{taxon}::{go}"] |= set(subs)
                sources["complexportal"]["entries"] += 1
                sources["complexportal"][f"taxon_{taxon}"] = 1

    if args.corum_file and os.path.exists(args.corum_file):
        for taxon, go, subs in parse_corum_file(args.corum_file):
            index[f"{taxon}::{go}"] |= set(subs)
            sources["corum"]["entries"] += 1
            sources["corum"][f"taxon_{taxon}"] = 1

    print("Source contribution:")
    for src, stats in sources.items():
        n_taxa = sum(1 for k in stats if k.startswith("taxon_"))
        print(f"  {src}: {stats['entries']} entries across {n_taxa} taxa")
    print(f"Total unique (taxon, complex GO term) keys: {len(index)}")

    out = {k: sorted(v) for k, v in index.items()}
    with open(args.output_file, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print(f"Wrote {args.output_file}")


if __name__ == "__main__":
    main()
