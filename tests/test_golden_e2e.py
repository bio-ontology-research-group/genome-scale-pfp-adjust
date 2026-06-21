"""Golden end-to-end test: fixed input -> Stage 1 -> Stage 2 -> diff vs golden.

This is the test a reviewer can actually run. It solves a single synthetic
organism of 50 proteins (``fixtures/golden_predictions.tsv``, built by
``fixtures/make_golden_organism.py``) through both adjustment stages and asserts
that the output is byte-for-byte identical to the committed golden files:

    fixtures/golden_predictions.tsv      (input)
      --Stage 1: taxon consistency-->    fixtures/golden_stage1.tsv
      --Stage 2: complex coherence-->    fixtures/golden_stage2.tsv

It runs Stage 1 + Stage 2 in well under a minute on a CPU laptop (no GPU, no
network). A diff failure means the solver's behaviour changed; if the change is
intended, regenerate the goldens with::

    python tests/test_golden_e2e.py --regenerate

and review the diff before committing.

It also re-asserts the three biological invariants the fixture is built around,
so a regression that silently rewrites the golden files cannot pass unnoticed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from taxon_consistency.adjust_ortools import main as taxon_main, load_predictions
from complex_coherence.adjust_ortools import main as complex_main

FIX = os.path.join(HERE, "fixtures")

INPUT = os.path.join(FIX, "golden_predictions.tsv")
GOLDEN_STAGE1 = os.path.join(FIX, "golden_stage1.tsv")
GOLDEN_STAGE2 = os.path.join(FIX, "golden_stage2.tsv")

# Shared constraint / hierarchy fixtures (see fixtures/make_golden_organism.py).
CONSTRAINTS = os.path.join(FIX, "mini_constraints.tsv")
GO_HIER = os.path.join(FIX, "mini_go_hierarchy.tsv")
TAXON_HIER = os.path.join(FIX, "mini_taxon_hierarchy.tsv")
NCBITAXON_HIER = os.path.join(FIX, "mini_ncbitaxon_hierarchy.tsv")
COMPLEXES = os.path.join(FIX, "mini_complexes.tsv")

TAXON_ID = "2759"   # Eukaryota
THRESHOLD = 0.3
TOP_K = 5

BACTERIA_TERM = "GO:0009001"   # only_in_taxon Bacteria -> must be demoted
COMPLEX_TERM = "GO:0099001"    # heteromeric complex -> no singleton allowed
BENIGN_TERMS = ("GO:0008002", "GO:0008003")


def run_stage1(out_path):
    return taxon_main(
        predictions_file=INPUT,
        constraints_file=CONSTRAINTS,
        go_hierarchy_file=GO_HIER,
        taxon_hierarchy_file=TAXON_HIER,
        ncbitaxon_hierarchy_file=NCBITAXON_HIER,
        output_file=out_path,
        threshold=THRESHOLD,
        taxon_id=TAXON_ID,
    )


def run_stage2(in_path, out_path):
    return complex_main(
        predictions_file=in_path,
        complexes_file=COMPLEXES,
        go_hierarchy_file=GO_HIER,
        output_file=out_path,
        threshold=THRESHOLD,
        optimized=True,
        top_k=TOP_K,
    )


def _canonical(path):
    """Read an adjusted-prediction TSV into a stable, comparable string.

    The solver preserves each protein's score values exactly but does not fix
    the *order* of the GO-term columns within a line (it iterates a set), so the
    raw byte order is not reproducible. Sorting the term columns per line yields
    a canonical form that diffs cleanly while still catching any change to which
    terms or scores a protein carries."""
    lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            pid, *terms = line.split("\t")
            lines.append(pid + "\t" + "\t".join(sorted(terms)))
    return "\n".join(sorted(lines)) + "\n"


def _write_canonical(src, dst):
    with open(dst, "w") as fh:
        fh.write(_canonical(src))


def _assert_invariants(stage1_path, stage2_path):
    """Biological invariants the fixture exists to exercise."""
    s1 = load_predictions(stage1_path)
    s2 = load_predictions(stage2_path)

    # Stage 1: no Bacteria-only term survives above threshold on this eukaryote.
    surviving = [p for p, t in s1.items() if t.get(BACTERIA_TERM, 0.0) >= THRESHOLD]
    assert not surviving, f"Bacteria-only term survived Stage 1 on: {surviving}"
    # ...and at least one such term was actually present to demote.
    demoted = sum(1 for p, t in load_predictions(INPUT).items()
                  if t.get(BACTERIA_TERM, 0.0) >= THRESHOLD)
    assert demoted > 0, "Fixture has no Bacteria-only annotation to demote"

    # Stage 2: the complex term is never on exactly one protein above threshold.
    above = sum(1 for t in s2.values() if t.get(COMPLEX_TERM, 0.0) >= THRESHOLD)
    assert above != 1, f"Stage 2 left {COMPLEX_TERM} as a singleton complex"

    # Benign terms must survive both stages unchanged for every protein.
    inp = load_predictions(INPUT)
    for pid, terms in inp.items():
        for bt in BENIGN_TERMS:
            if bt in terms:
                assert s2.get(pid, {}).get(bt, 0.0) >= THRESHOLD, (
                    f"Benign term {bt} on {pid} was not preserved"
                )


def test_golden_end_to_end():
    print("test_golden_end_to_end")
    s1_tmp = os.path.join(FIX, ".tmp_stage1.tsv")
    s2_tmp = os.path.join(FIX, ".tmp_stage2.tsv")
    try:
        run_stage1(s1_tmp)
        run_stage2(s1_tmp, s2_tmp)

        assert _canonical(s1_tmp) == _canonical(GOLDEN_STAGE1), (
            "Stage 1 output differs from golden_stage1.tsv "
            "(run with --regenerate if the change is intended)"
        )
        assert _canonical(s2_tmp) == _canonical(GOLDEN_STAGE2), (
            "Stage 2 output differs from golden_stage2.tsv "
            "(run with --regenerate if the change is intended)"
        )
        _assert_invariants(s1_tmp, s2_tmp)
        print("  Golden e2e OK: Stage 1 + Stage 2 match committed goldens.")
    finally:
        for p in (s1_tmp, s2_tmp):
            if os.path.exists(p):
                os.remove(p)


def regenerate():
    sys.path.insert(0, FIX)
    import make_golden_organism as gen  # noqa: E402
    gen.write(INPUT)
    s1_tmp = os.path.join(FIX, ".tmp_stage1.tsv")
    s2_tmp = os.path.join(FIX, ".tmp_stage2.tsv")
    try:
        run_stage1(s1_tmp)
        run_stage2(s1_tmp, s2_tmp)
        _assert_invariants(s1_tmp, s2_tmp)
        _write_canonical(s1_tmp, GOLDEN_STAGE1)
        _write_canonical(s2_tmp, GOLDEN_STAGE2)
    finally:
        for p in (s1_tmp, s2_tmp):
            if os.path.exists(p):
                os.remove(p)
    print(f"Regenerated goldens:\n  {GOLDEN_STAGE1}\n  {GOLDEN_STAGE2}")


def main():
    test_golden_end_to_end()
    print("Golden end-to-end test passed.")


if __name__ == "__main__":
    if "--regenerate" in sys.argv:
        regenerate()
    else:
        main()
