"""
Validate Stage-2 (complex coherence) promotions/demotions against curated
ComplexPortal/CORUM subunit lists.

Stage 2 repairs each singleton-complex violation by either DEMOTING the term
from the lone protein or PROMOTING it onto an additional protein, chosen by
the solver on cost (score margin), not biology. This script closes that loop
by diffing the Stage-1-only output against the Stage-1+2 output (the diff is,
by construction, confined to complex terms and their GO ancestors) and scoring
each promotion/demotion against a curated subunit index:

  promotion precision (covered) = fraction of promoted (taxon, complex, protein)
      triples, restricted to those with a curated reference, where the promoted
      protein IS a curated subunit of the complex in that organism.
  demotion good-rate (covered)  = fraction of demoted singletons, restricted to
      those with a curated reference, where the demoted protein is NOT a curated
      subunit (i.e. the lone prediction was an unsupported false positive, so
      demotion was the right call).

Coverage is reported alongside every precision: a missing index key means "no
curated reference for this (organism, complex)", not "no subunits exist".

Subunit index format (from build_subunit_index.py): JSON dict keyed
"<taxon_id>::<GO term>" -> sorted list of curated UniProt accessions.
"""

import argparse
import glob
import json
import os
import re
from collections import defaultdict

FNAME_RE = re.compile(r"optimized_fold_(\d+)_taxon_(\d+)\.tsv$")


def protein_to_accession(pid):
    """`A0A6P5QAZ4_MUSCR` -> `A0A6P5QAZ4`; a bare accession is returned as-is."""
    return pid.split("_", 1)[0] if "_" in pid else pid


def load_terms_above(path, threshold):
    """protein_id -> set(GO terms with score > threshold)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            prot = parts[0]
            terms = set()
            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go, score = tok.split("|", 1)
                try:
                    if float(score) > threshold:
                        terms.add(go)
                except ValueError:
                    continue
            out[prot] = terms
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage1_dir", required=True,
                    help="Stage-1-only optimized dir (input to Stage 2)")
    ap.add_argument("--stage2_dir", required=True,
                    help="Stage-1+2 (complex coherence) optimized dir")
    ap.add_argument("--subunit_index", required=True, help="subunit_index.json")
    ap.add_argument("--threshold", type=float, required=True,
                    help="tau used for this predictor (DeepGO-SE 0.3, else 0.001)")
    ap.add_argument("--label", default="", help="predictor label for the report")
    ap.add_argument("--output_file", default=None, help="optional JSON report path")
    args = ap.parse_args()

    with open(args.subunit_index) as fh:
        index = json.load(fh)
    index = {k: set(v) for k, v in index.items()}

    promo_total = promo_cov = promo_tp = 0
    demo_total = demo_cov = demo_good = 0
    promo_examples, demo_harmful_examples = [], []
    per_taxon = defaultdict(lambda: {"promo": 0, "promo_cov": 0, "promo_tp": 0})

    s2_files = sorted(glob.glob(os.path.join(args.stage2_dir, "optimized_fold_*_taxon_*.tsv")))
    n_pairs = 0
    for s2 in s2_files:
        m = FNAME_RE.search(os.path.basename(s2))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        s1 = os.path.join(args.stage1_dir, f"optimized_fold_{fold}_taxon_{taxon}.tsv")
        if not os.path.exists(s1):
            continue
        n_pairs += 1
        t1 = load_terms_above(s1, args.threshold)
        t2 = load_terms_above(s2, args.threshold)
        for prot in set(t1) | set(t2):
            a1 = t1.get(prot, set())
            a2 = t2.get(prot, set())
            acc = protein_to_accession(prot)
            for go in (a2 - a1):  # promotions
                promo_total += 1
                key = f"{taxon}::{go}"
                per_taxon[taxon]["promo"] += 1
                if key in index:
                    promo_cov += 1
                    per_taxon[taxon]["promo_cov"] += 1
                    if acc in index[key]:
                        promo_tp += 1
                        per_taxon[taxon]["promo_tp"] += 1
                        if len(promo_examples) < 25:
                            promo_examples.append([taxon, go, acc, "subunit"])
            for go in (a1 - a2):  # demotions
                demo_total += 1
                key = f"{taxon}::{go}"
                if key in index:
                    demo_cov += 1
                    if acc not in index[key]:
                        demo_good += 1
                    elif len(demo_harmful_examples) < 25:
                        demo_harmful_examples.append([taxon, go, acc, "curated_subunit_removed"])

    def pct(a, b):
        return round(100.0 * a / b, 2) if b else None

    report = {
        "label": args.label,
        "n_organism_pairs": n_pairs,
        "threshold": args.threshold,
        "promotions": {
            "total": promo_total,
            "with_curated_reference": promo_cov,
            "coverage_pct": pct(promo_cov, promo_total),
            "precision_pct_on_covered": pct(promo_tp, promo_cov),
            "true_subunit_promotions": promo_tp,
        },
        "demotions": {
            "total": demo_total,
            "with_curated_reference": demo_cov,
            "coverage_pct": pct(demo_cov, demo_total),
            "good_demotion_pct_on_covered": pct(demo_good, demo_cov),
            "harmful_demotions_curated_subunit_removed": demo_cov - demo_good,
        },
        "promotion_examples": promo_examples,
        "harmful_demotion_examples": demo_harmful_examples,
    }
    print(json.dumps(report, indent=2))
    if args.output_file:
        with open(args.output_file, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"\nWrote {args.output_file}")


if __name__ == "__main__":
    main()
