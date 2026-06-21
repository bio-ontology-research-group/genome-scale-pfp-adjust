#!/usr/bin/env python3
"""
Materialize the obligate heteromeric complex set used by Stage 2.

The set is the descendants of GO:0032991 (protein-containing complex) in the
cellular-component hierarchy, minus the hand-classified homodimer/auxiliary
terms (classification 'h' or 'a' in complex_coherence/protein_complexes.tsv).
This is the same computation performed at runtime by
complex_coherence/adjust_ortools.py::get_heteromeric_complexes; here we write it
to a static, checksummed artefact for reproducibility (paper Section 2.1).

Usage:
    python data_prep/materialize_complex_list.py \
        --go_hierarchy_cc data/go_hierarchy_cc.tsv \
        --complexes complex_coherence/protein_complexes.tsv \
        --output data/constraints/heteromeric_complexes_2025_10.tsv

Inputs are the GO release 2025_10 CC hierarchy (child<TAB>parent) and the
curated classification file. Output is one GO term per line, sorted, with a
short provenance header. Run `sha256sum` on the output to obtain the checksum
cited in the manuscript.
"""

import argparse
import hashlib
from collections import defaultdict

ROOT = "GO:0032991"  # protein-containing complex
GO_RELEASE = "2025_10"


def load_child_to_parents(path):
    """Load child<TAB>parent rows into {child: {parents}} (load_go_hierarchy)."""
    hier = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and parts[1]:
                hier[parts[0]].add(parts[1])
    return hier


def all_descendants(root, child_to_parents):
    """Every term reachable below root, matching get_all_children's BFS."""
    children_of = defaultdict(set)
    for child, parents in child_to_parents.items():
        for parent in parents:
            children_of[parent].add(child)
    seen, queue = set(), [root]
    while queue:
        cur = queue.pop(0)
        for child in sorted(children_of.get(cur, ())):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return seen


def load_homodimer_terms(path):
    """GO terms classified homodimer ('h') or auxiliary/ambiguous ('a')."""
    homo = set()
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            row = line.rstrip("\n").split("\t")
            if len(row) >= 2 and row[1] in ("h", "a"):
                homo.add(row[0])
    return homo


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--go_hierarchy_cc", default="data/go_hierarchy_cc.tsv")
    ap.add_argument("--complexes", default="complex_coherence/protein_complexes.tsv")
    ap.add_argument("--output", default="data/constraints/heteromeric_complexes_2025_10.tsv")
    args = ap.parse_args()

    hier = load_child_to_parents(args.go_hierarchy_cc)
    descendants = all_descendants(ROOT, hier)
    descendants.add(ROOT)
    homo = load_homodimer_terms(args.complexes)
    heteromeric = sorted(descendants - homo)

    with open(args.output, "w") as fh:
        fh.write(f"# Obligate heteromeric complex terms (GAEF set), GO release {GO_RELEASE}\n")
        fh.write(f"# = descendants({ROOT}) in the CC hierarchy minus {len(homo)} "
                 f"hand-classified homodimer/auxiliary terms\n")
        fh.write(f"# count\t{len(heteromeric)}\n")
        for go in heteromeric:
            fh.write(go + "\n")

    with open(args.output, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()

    print(f"descendants({ROOT})+self: {len(descendants)}")
    print(f"homodimer/auxiliary terms: {len(homo)}")
    print(f"heteromeric obligate complexes: {len(heteromeric)}")
    print(f"wrote {args.output}")
    print(f"sha256: {digest}")


if __name__ == "__main__":
    main()
