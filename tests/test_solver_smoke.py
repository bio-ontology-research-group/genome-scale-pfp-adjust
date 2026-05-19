"""
End-to-end smoke test for both solver stages.

Stage 1 (taxon consistency):
    Predictions contain a bacteria-only GO term (GO:0009001) on a protein in a
    eukaryotic organism (taxon NCBITaxon_2759). The solver must demote it.

Stage 2 (complex coherence):
    Predictions contain a heteromeric complex GO term (GO:0099001) above
    threshold on exactly one protein. The solver must repair this either by
    demoting from the singleton or by promoting a second protein.

Both stages should finish in under 5 s on a laptop.
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from taxon_consistency.adjust_ortools import (
    main as taxon_main,
    load_predictions,
)
from complex_coherence.adjust_ortools import main as complex_main


FIX = os.path.join(HERE, "fixtures")


def test_stage1_taxon_consistency():
    print("test_stage1_taxon_consistency")
    with tempfile.NamedTemporaryFile(suffix="_taxon.tsv", delete=False) as f:
        out_path = f.name
    try:
        total_flips, num_annotated = taxon_main(
            predictions_file=os.path.join(FIX, "mini_predictions_taxon.tsv"),
            constraints_file=os.path.join(FIX, "mini_constraints.tsv"),
            go_hierarchy_file=os.path.join(FIX, "mini_go_hierarchy.tsv"),
            taxon_hierarchy_file=os.path.join(FIX, "mini_taxon_hierarchy.tsv"),
            ncbitaxon_hierarchy_file=os.path.join(FIX, "mini_ncbitaxon_hierarchy.tsv"),
            output_file=out_path,
            threshold=0.3,
            taxon_id="2759",
        )
        assert num_annotated > 0, "Expected some annotated predictions"
        adjusted = load_predictions(out_path)

        protein1 = adjusted.get("protein1", {})
        # GO:0009001 was bacteria-only on a eukaryote: solver should have demoted it
        if "GO:0009001" in protein1:
            assert protein1["GO:0009001"] < 0.3, (
                f"Expected GO:0009001 to be demoted below 0.3 on protein1, got {protein1['GO:0009001']}"
            )
        # Non-violating annotations should be preserved
        assert adjusted.get("protein2", {}).get("GO:0008002", 0.0) >= 0.3, \
            "Expected GO:0008002 on protein2 to be preserved"
        print(f"  Stage 1 OK: total_flips={total_flips}, num_annotated={num_annotated}")
    finally:
        if os.path.exists(out_path):
            os.remove(out_path)


def test_stage2_complex_coherence():
    print("test_stage2_complex_coherence")
    with tempfile.NamedTemporaryFile(suffix="_complex.tsv", delete=False) as f:
        out_path = f.name
    try:
        total_flips, num_annotated = complex_main(
            predictions_file=os.path.join(FIX, "mini_predictions_complex.tsv"),
            complexes_file=os.path.join(FIX, "mini_complexes.tsv"),
            go_hierarchy_file=os.path.join(FIX, "mini_go_hierarchy.tsv"),
            output_file=out_path,
            threshold=0.3,
            optimized=True,
            top_k=5,
        )
        adjusted = load_predictions(out_path)
        p1 = adjusted.get("protein1", {}).get("GO:0099001", 0.0)
        p2 = adjusted.get("protein2", {}).get("GO:0099001", 0.0)
        p3 = adjusted.get("protein3", {}).get("GO:0099001", 0.0)
        # After repair, GO:0099001 must NOT be on exactly one protein above threshold.
        above = sum(1 for x in (p1, p2, p3) if x >= 0.3)
        assert above != 1, (
            f"Expected GO:0099001 to be on 0 or >=2 proteins above threshold, "
            f"got {above} (p1={p1}, p2={p2}, p3={p3})"
        )
        print(f"  Stage 2 OK: total_flips={total_flips}, num_annotated={num_annotated}")
    finally:
        if os.path.exists(out_path):
            os.remove(out_path)


def main():
    test_stage1_taxon_consistency()
    test_stage2_complex_coherence()
    print("All solver smoke tests passed.")


if __name__ == "__main__":
    main()
