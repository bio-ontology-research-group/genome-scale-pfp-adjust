"""Deterministically (re)generate the golden end-to-end input fixture.

Writes ``golden_predictions.tsv``: a synthetic eukaryotic organism
(NCBITaxon_2759) of 50 proteins used by ``tests/test_golden_e2e.py``. The
file is fully deterministic (no RNG), so re-running this script must not change
its contents. The constraint/hierarchy fixtures it is solved against are the
shared ``mini_*.tsv`` files in this directory:

  * GO:0009001 is ``only_in_taxon NCBITaxon_2`` (Bacteria) and is disjoint from
    the organism's lineage, so every protein carrying it above threshold must be
    demoted by Stage 1 (taxon consistency).
  * GO:0099001 is an obligate heteromeric complex term; the input places it
    above threshold on exactly one protein (a singleton), which Stage 2
    (complex coherence) must repair.
  * GO:0008002 / GO:0008003 are benign terms with no taxon or complex
    constraint and must survive both stages unchanged.

To regenerate the input and the golden outputs together, run the helper at the
bottom of ``tests/test_golden_e2e.py`` (``python tests/test_golden_e2e.py
--regenerate``) rather than this script alone.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

N_PROTEINS = 50
# Every BACTERIA_EVERY-th protein (1-indexed) carries the Bacteria-only term.
BACTERIA_EVERY = 4


def _benign_scores(i):
    """Two benign terms, both deterministically above the 0.3 threshold."""
    s2 = 0.50 + ((i * 7) % 5) * 0.08   # in {0.50, 0.58, 0.66, 0.74, 0.82}
    s3 = 0.45 + ((i * 3) % 4) * 0.10   # in {0.45, 0.55, 0.65, 0.75}
    return s2, s3


def build_rows():
    rows = []
    for i in range(1, N_PROTEINS + 1):
        pid = f"p{i:02d}"
        s2, s3 = _benign_scores(i)
        terms = [f"GO:0008002|{s2:.2f}", f"GO:0008003|{s3:.2f}"]
        # Bacteria-only term on a eukaryote -> Stage 1 must demote it.
        if i % BACTERIA_EVERY == 1:
            terms.insert(0, "GO:0009001|0.90")
        # Heteromeric complex term: above threshold on exactly one protein
        # (p01) -> Stage 2 must repair the singleton. p02 is the cheapest
        # below-threshold promotion candidate; p03 is a more expensive one.
        if i == 1:
            terms.append("GO:0099001|0.85")
        elif i == 2:
            terms.append("GO:0099001|0.28")
        elif i == 3:
            terms.append("GO:0099001|0.12")
        rows.append((pid, terms))
    return rows


def write(path=None):
    path = path or os.path.join(HERE, "golden_predictions.tsv")
    rows = build_rows()
    with open(path, "w") as fh:
        for pid, terms in rows:
            fh.write(pid + "\t" + "\t".join(terms) + "\n")
    return path


if __name__ == "__main__":
    p = write()
    print(f"Wrote {p} ({N_PROTEINS} proteins)")
