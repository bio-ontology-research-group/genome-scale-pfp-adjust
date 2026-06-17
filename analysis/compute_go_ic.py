"""
Precompute GO-term information content (IC) from the ground-truth proteome
annotations, for the IC-weighted cost in the cost-function ablation.

IC(t) = -log2( freq(t) / N ), where freq(t) is the number of proteins annotated
with t (after propagating each protein's terms to all GO ancestors, true-path
rule) and N is the number of annotated proteins. Writes go_ic.tsv (GO<TAB>IC).
"""
import argparse
import csv
import glob
import math
import os
from collections import defaultdict


def load_go_parents(path):
    parents = defaultdict(set)
    with open(path) as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) >= 2:
                parents[row[0]].add(row[1])
    return parents


def closure(seed, parents):
    out, stack = set(seed), list(seed)
    while stack:
        for p in parents.get(stack.pop(), ()):
            if p not in out:
                out.add(p)
                stack.append(p)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotations_dir", required=True)
    ap.add_argument("--go_hierarchy_file", default="data/go_hierarchy.tsv")
    ap.add_argument("--output_file", required=True)
    args = ap.parse_args()

    parents = load_go_parents(args.go_hierarchy_file)
    cnt = defaultdict(int)
    n_prot = 0
    for f in glob.glob(os.path.join(args.annotations_dir, "annots_taxon_*.tsv")):
        for line in open(f):
            parts = line.rstrip("\n").split("\t")
            terms = [t for t in parts[1:] if t.startswith("GO:")]
            if not terms:
                continue
            n_prot += 1
            for t in closure(set(terms), parents):
                cnt[t] += 1
    with open(args.output_file, "w") as out:
        for t, c in cnt.items():
            ic = -math.log2(c / n_prot) if c > 0 else 0.0
            out.write(f"{t}\t{ic:.6f}\n")
    print(f"Wrote IC for {len(cnt)} terms over {n_prot} annotated proteins to {args.output_file}")


if __name__ == "__main__":
    main()
