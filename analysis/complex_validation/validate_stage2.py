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
  Demotions are separated into listed and unlisted proteins. Absence from a
      curated list is not a validated negative or evidence of a correct removal.

Coverage is reported alongside every precision: a missing index key means "no
curated reference for this (organism, complex)", not "no subunits exist".

Subunit index format (from build_subunit_index.py): JSON dict keyed
"<taxon_id>::<GO term>" -> sorted list of curated UniProt accessions.
"""

import argparse
import glob
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gaef_patches.ontology import Ontology
from complex_coherence.adjust_ortools import load_homodimer_terms, get_heteromeric_complexes

FNAME_RE = re.compile(r"optimized_fold_(\d+)_taxon_(\d+)\.tsv$")


def file_sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def protein_to_accession(pid, mapping=None):
    """Resolve entry names explicitly; a gene mnemonic is not an accession."""
    if mapping and pid in mapping:
        return mapping[pid]
    candidate = pid.split('_', 1)[0]
    if re.fullmatch(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})', candidate):
        return candidate
    return None


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
    ap.add_argument("--protein-id-map", help="JSON mapping UniProt entry names to accessions")
    ap.add_argument("--go", required=True, help="Exact GO OBO used for the complex repair")
    ap.add_argument("--complex-exclusions", default=str(ROOT/'complex_coherence/protein_complexes.tsv'))
    args = ap.parse_args()

    with open(args.subunit_index) as fh:
        index = json.load(fh)
    if not index:
        raise ValueError('Empty curated reference index')
    go = Ontology(args.go, with_rels=False)
    complex_terms = get_heteromeric_complexes(go.normalize_terms(load_homodimer_terms(args.complex_exclusions)), set(),
        {term: set(obj['is_a']) for term, obj in go.ont.items()})
    raw_reference_count = len(index)
    # Complex Portal also supplies process/function/location annotations. Those
    # references do not establish membership of a complex class.
    normalized_index = defaultdict(set)
    for key, accessions in index.items():
        taxon, term = key.split('::', 1)
        term = go.canonical_id(term)
        if term in complex_terms:
            normalized_index[taxon+'::'+term].update(accessions)
    index = normalized_index
    mapping = json.load(open(args.protein_id_map)) if args.protein_id_map else {}
    unmapped_promotions = unmapped_demotions = 0

    promo_total = promo_cov = promo_tp = 0
    demo_total = demo_cov = demo_good = 0
    promo_examples, demo_harmful_examples = [], []
    per_taxon = defaultdict(lambda: {"promo": 0, "promo_cov": 0, "promo_tp": 0})

    s2_files = sorted(glob.glob(os.path.join(args.stage2_dir, "optimized_fold_*_taxon_*.tsv")))
    if not s2_files:
        raise FileNotFoundError('No Stage 2 files found')
    s1_names = {os.path.basename(p) for p in glob.glob(os.path.join(args.stage1_dir, 'optimized_fold_*_taxon_*.tsv'))}
    if s1_names != {os.path.basename(p) for p in s2_files}:
        raise ValueError('Stage 1 and Stage 2 organism manifests differ')
    n_pairs = 0
    pair_inputs = []
    for s2 in s2_files:
        m = FNAME_RE.search(os.path.basename(s2))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        s1 = os.path.join(args.stage1_dir, f"optimized_fold_{fold}_taxon_{taxon}.tsv")
        if not os.path.exists(s1):
            continue
        pair_inputs.append({"taxon": taxon, "stage1": file_sha(s1), "stage2": file_sha(s2)})
        n_pairs += 1
        t1 = load_terms_above(s1, args.threshold)
        t2 = load_terms_above(s2, args.threshold)
        for prot in sorted(set(t1) | set(t2)):
            a1 = t1.get(prot, set())
            a2 = t2.get(prot, set())
            acc = protein_to_accession(prot, mapping)
            for go in sorted(a2 - a1):  # promotions
                promo_total += 1
                key = f"{taxon}::{go}"
                per_taxon[taxon]["promo"] += 1
                if key in index and acc is None:
                    unmapped_promotions += 1
                if key in index and acc is not None:
                    promo_cov += 1
                    per_taxon[taxon]["promo_cov"] += 1
                    if acc in index[key]:
                        promo_tp += 1
                        per_taxon[taxon]["promo_tp"] += 1
                        if len(promo_examples) < 25:
                            promo_examples.append([taxon, go, acc, "subunit"])
            for go in sorted(a1 - a2):  # demotions
                demo_total += 1
                key = f"{taxon}::{go}"
                if key in index and acc is None:
                    unmapped_demotions += 1
                if key in index and acc is not None:
                    demo_cov += 1
                    if acc not in index[key]:
                        demo_good += 1
                    elif len(demo_harmful_examples) < 25:
                        demo_harmful_examples.append([taxon, go, acc, "curated_subunit_removed"])

    def pct(a, b):
        return round(100.0 * a / b, 2) if b else None

    report = {
        "provenance": {"script_sha256": file_sha(__file__), "index_sha256": file_sha(args.subunit_index),
            "go_sha256": file_sha(args.go), "complex_exclusions_sha256": file_sha(args.complex_exclusions),
            "raw_reference_keys": raw_reference_count, "complex_reference_keys": len(index),
            "mapping_sha256": file_sha(args.protein_id_map) if args.protein_id_map else None, "pairs": pair_inputs},
        "label": args.label,
        "n_organism_pairs": n_pairs,
        "threshold": args.threshold,
        "interpretation": "Absence from a curated membership list is not proof of non-membership.",
        "unmapped_promotions_with_reference": unmapped_promotions,
        "unmapped_demotions_with_reference": unmapped_demotions,
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
            "unlisted_demotion_pct_on_covered": pct(demo_good, demo_cov),
            "curated_subunit_demotions": demo_cov - demo_good,
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
