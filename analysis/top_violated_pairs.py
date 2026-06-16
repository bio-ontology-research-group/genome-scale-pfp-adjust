"""
Scan the curated ground-truth proteome annotations for taxon-constraint
violations and report the most frequently violated (constrained GO term,
constraint) pairs across organisms --- the per-pair breakdown behind the
aggregate 6.2% (32/514) figure.

A curated annotation violates a constraint when an annotated term's GO closure
(the term or any ancestor) contains a constrained term h such that:
  - never_in_taxon(x) with x in the organism's lineage, or
  - only_in_taxon(...) none of which is in the organism's lineage.

Reuses the solver's loaders so semantics match the CP-SAT / naive-filter paths.
"""
import argparse
import glob
import os
import re
from collections import defaultdict

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from taxon_consistency.adjust_ortools import (  # noqa: E402
    load_constraints, load_go_hierarchy, load_taxon_hierarchy,
)

FN = re.compile(r"annots_taxon_(\d+)\.tsv$")


def transitive(seed, adj):
    out, stack = set(seed), list(seed)
    while stack:
        for nxt in adj.get(stack.pop(), ()):
            if nxt not in out:
                out.add(nxt)
                stack.append(nxt)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations_dir", required=True)
    ap.add_argument("--constraints_file", default="data/go_taxon_constraints_extracted_obo.tsv")
    ap.add_argument("--go_hierarchy_file", default="data/go_hierarchy.tsv")
    ap.add_argument("--taxon_hierarchy_file", default="data/taxon_hierarchy.tsv")
    ap.add_argument("--ncbitaxon_hierarchy_file", default="data/ncbitaxon_hierarchy.tsv")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--output_file", default=None)
    args = ap.parse_args()

    in_taxon, never_in = load_constraints(args.constraints_file)
    go_parents = load_go_hierarchy(args.go_hierarchy_file)
    tax_isa, _, _ = load_taxon_hierarchy(args.taxon_hierarchy_file, args.ncbitaxon_hierarchy_file)

    pair_orgs = defaultdict(set)      # (h, type, ctaxa) -> set(org taxon ids)
    pair_annots = defaultdict(int)    # (h, type, ctaxa) -> annotation count
    n_orgs = n_incons = 0

    for f in sorted(glob.glob(os.path.join(args.annotations_dir, "annots_taxon_*.tsv"))):
        m = FN.search(os.path.basename(f))
        if not m:
            continue
        org = m.group(1)
        taxon = f"NCBITaxon_{org}"
        lineage = transitive({taxon}, tax_isa)
        # constrained terms violated for this organism
        violated = {}
        for g, taxa in never_in.items():
            hit = [t for t in taxa if t in lineage]
            if hit:
                violated[g] = ("never_in_taxon", ",".join(sorted(hit)))
        for g, taxa in in_taxon.items():
            if taxa and not any(t in lineage for t in taxa):
                violated[g] = ("only_in_taxon", ",".join(sorted(taxa)))
        if not violated:
            n_orgs += 1
            continue
        org_has_violation = False
        for line in open(f):
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            terms = [t for t in parts[1:] if t.startswith("GO:")]
            if not terms:
                continue
            closure = transitive(set(terms), go_parents)
            for h in (closure & set(violated)):
                typ, ctaxa = violated[h]
                key = (h, typ, ctaxa)
                pair_orgs[key].add(org)
                pair_annots[key] += 1
                org_has_violation = True
        n_orgs += 1
        if org_has_violation:
            n_incons += 1

    ranked = sorted(pair_orgs, key=lambda k: (len(pair_orgs[k]), pair_annots[k]), reverse=True)
    print(f"# organisms scanned: {n_orgs}; with >=1 curated violation: {n_incons} "
          f"({100.0*n_incons/n_orgs:.2f}%)")
    print(f"# distinct violated (GO term, constraint, taxon) pairs: {len(ranked)}")
    print("\nrank\tGO_term\tconstraint\tconstraint_taxa\tn_organisms\tn_annotations")
    lines = []
    for i, k in enumerate(ranked[:args.top], 1):
        h, typ, ctaxa = k
        row = f"{i}\t{h}\t{typ}\t{ctaxa}\t{len(pair_orgs[k])}\t{pair_annots[k]}"
        print(row)
        lines.append(row)
    if args.output_file:
        with open(args.output_file, "w") as fh:
            fh.write(f"# organisms_scanned={n_orgs} inconsistent={n_incons}\n")
            fh.write("rank\tGO_term\tconstraint\tconstraint_taxa\tn_organisms\tn_annotations\n")
            fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
