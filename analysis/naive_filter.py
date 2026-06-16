"""
Naive taxon filter baseline (the "why not a SQL query?" comparison).

Given the organism's taxon, drop predicted GO terms that violate its taxon
constraints --- no CP-SAT, no joint inference, no disjoint/union taxon algebra,
no complex coherence. Two variants:

  --no-propagate : drop a predicted term only if the term ITSELF carries a
                   violated taxon constraint.
  (default)      : also drop predicted terms whose GO closure (the term or any
                   ancestor) is a violated constrained term --- i.e. propagate
                   violations down the GO hierarchy (true-path rule).

Output TSVs keep the original predictions filename so they can be fed straight
into pipeline/evaluate_directory.py, whose own consistency checker
(`satisfiable`) is then the independent judge of whether the filter achieved
genome-scale taxon consistency.

Reuses the solver's loaders (taxon_consistency.adjust_ortools) so constraint /
hierarchy / ID parsing is identical to the CP-SAT path.
"""
import argparse
import glob
import os
import re
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from taxon_consistency.adjust_ortools import (  # noqa: E402
    load_constraints, load_go_hierarchy, load_taxon_hierarchy,
)

FN = re.compile(r"_fold_(\d+)_taxon_(\d+)\.tsv$")


def transitive(seed, adj):
    """Transitive closure of `seed` under adjacency map `adj` (node -> set)."""
    out, stack = set(seed), list(seed)
    while stack:
        for nxt in adj.get(stack.pop(), ()):
            if nxt not in out:
                out.add(nxt)
                stack.append(nxt)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions_dir", required=True)
    ap.add_argument("--output_dir", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--constraints_file", default="data/go_taxon_constraints_extracted_obo.tsv")
    ap.add_argument("--go_hierarchy_file", default="data/go_hierarchy.tsv")
    ap.add_argument("--taxon_hierarchy_file", default="data/taxon_hierarchy.tsv")
    ap.add_argument("--ncbitaxon_hierarchy_file", default="data/ncbitaxon_hierarchy.tsv")
    ap.add_argument("--no-propagate", dest="propagate", action="store_false")
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    in_taxon, never_in = load_constraints(args.constraints_file)
    go_parents = load_go_hierarchy(args.go_hierarchy_file)          # child -> {parents}
    go_children = {}                                                 # parent -> {children}
    for c, ps in go_parents.items():
        for p in ps:
            go_children.setdefault(p, set()).add(c)
    tax_isa, _, _ = load_taxon_hierarchy(args.taxon_hierarchy_file, args.ncbitaxon_hierarchy_file)

    for pred in sorted(glob.glob(os.path.join(args.predictions_dir, "predictions_fold_*_taxon_*.tsv"))):
        m = FN.search(os.path.basename(pred))
        if not m:
            continue
        taxon = f"NCBITaxon_{m.group(2)}"
        lineage = transitive({taxon}, tax_isa)        # org taxon + all ancestors

        # constrained terms violated for this lineage
        violated = set()
        for g, taxa in never_in.items():
            if any(t in lineage for t in taxa):
                violated.add(g)
        for g, taxa in in_taxon.items():
            if taxa and not any(t in lineage for t in taxa):
                violated.add(g)
        # propagate: a violated term forbids all its GO descendants (true-path)
        forbidden = transitive(violated, go_children) if args.propagate else set(violated)

        out_path = os.path.join(args.output_dir, os.path.basename(pred))
        with open(pred) as fin, open(out_path, "w") as fout:
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0]:
                    continue
                kept = [parts[0]]
                for tok in parts[1:]:
                    if "|" not in tok:
                        continue
                    go, sc = tok.split("|", 1)
                    try:
                        above = float(sc) > args.threshold
                    except ValueError:
                        above = False
                    # drop above-threshold predictions that are forbidden; keep the rest as-is
                    if above and go in forbidden:
                        continue
                    kept.append(tok)
                fout.write("\t".join(kept) + "\n")
    print(f"Naive filter (propagate={args.propagate}) wrote to {args.output_dir}")


if __name__ == "__main__":
    main()
