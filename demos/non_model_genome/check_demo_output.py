#!/usr/bin/env python3
"""Check the non-model genome demo output.

The demo uses a small prediction table for Thermococcus kodakarensis KOD1
(NCBITaxon_69014). It deliberately includes two archaeon-incompatible
biological-process terms and one singleton heteromeric-complex term so both
solver stages have visible work to do.
"""

from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "data" / "predictions" / "predictions_fold_00_taxon_69014.tsv"
OUTPUT = ROOT / "output" / "optimized_fold_00_taxon_69014.tsv"
THRESHOLD = 0.3

TAXON_INCOMPATIBLE = {
    "GO:0032501": "multicellular organismal process",
    "GO:0045087": "innate immune response",
}
COMPLEX_TERM = "GO:0009333"


def load_predictions(path):
    predictions = {}
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            protein, *term_scores = line.split("\t")
            predictions[protein] = {}
            for term_score in term_scores:
                term, score = term_score.split("|")
                predictions[protein][term] = float(score)
    return predictions


def count_above(predictions):
    return sum(
        1
        for terms in predictions.values()
        for score in terms.values()
        if score > THRESHOLD
    )


def proteins_above(predictions, term):
    return sorted(
        protein
        for protein, terms in predictions.items()
        if terms.get(term, 0.0) > THRESHOLD
    )


def main():
    if not OUTPUT.exists():
        raise SystemExit(f"Missing demo output: {OUTPUT}")

    original = load_predictions(INPUT)
    adjusted = load_predictions(OUTPUT)

    removed = []
    for term, label in TAXON_INCOMPATIBLE.items():
        before = proteins_above(original, term)
        after = proteins_above(adjusted, term)
        if before and not after:
            removed.append(f"{term} ({label}) from {', '.join(before)}")
        if after:
            raise AssertionError(f"{term} still above threshold after adjustment: {after}")

    complex_before = proteins_above(original, COMPLEX_TERM)
    complex_after = proteins_above(adjusted, COMPLEX_TERM)
    if len(complex_before) != 1:
        raise AssertionError(f"Expected singleton complex before adjustment, got {complex_before}")
    if len(complex_after) == 1:
        raise AssertionError(f"Singleton complex survived adjustment: {complex_after}")

    print("Demo genome: Thermococcus kodakarensis KOD1 (NCBITaxon_69014)")
    print(f"Proteins: {len(original)}")
    print(f"Above-threshold annotations before: {count_above(original)}")
    print(f"Above-threshold annotations after: {count_above(adjusted)}")
    print("Taxon-stage removals:")
    for item in sorted(removed):
        print(f"  - {item}")
    print(f"Complex term {COMPLEX_TERM} before: {', '.join(complex_before)}")
    print(f"Complex term {COMPLEX_TERM} after: {', '.join(complex_after)}")
    print("OK: taxon consistency and complex coherence invariants hold.")


if __name__ == "__main__":
    main()
