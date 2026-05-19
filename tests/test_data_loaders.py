"""
Smoke test for the bundled data loaders.

Verifies that each TSV in data/ parses with the loader functions and produces
well-formed, non-empty dicts.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from taxon_consistency.adjust_ortools import (
    load_constraints,
    load_go_hierarchy,
    load_taxon_hierarchy,
)


def test_load_constraints():
    path = os.path.join(ROOT, "data", "go_taxon_constraints_updated.tsv")
    in_tax, never_in_tax = load_constraints(path)
    assert len(in_tax) > 0, f"No only_in_taxon constraints in {path}"
    assert len(never_in_tax) > 0, f"No never_in_taxon constraints in {path}"
    sample = next(iter(in_tax))
    assert sample.startswith("GO:"), f"Expected GO: prefix, got {sample!r}"
    print(f"  load_constraints: {len(in_tax)} only_in_taxon, {len(never_in_tax)} never_in_taxon")


def test_load_go_hierarchy():
    path = os.path.join(ROOT, "data", "go_hierarchy.tsv")
    hier = load_go_hierarchy(path)
    assert len(hier) > 1000, f"GO hierarchy looks truncated: {len(hier)} edges"
    sample_key = next(iter(hier))
    assert sample_key.startswith("GO:"), f"Expected GO: prefix, got {sample_key!r}"
    print(f"  load_go_hierarchy: {len(hier)} GO children")


def test_load_taxon_hierarchy():
    path = os.path.join(ROOT, "data", "taxon_hierarchy.tsv")
    is_a, disjoint, union = load_taxon_hierarchy(path, ncbitaxon_hierarchy_file=None)
    assert len(is_a) > 0 or len(union) > 0, "Empty taxon hierarchy"
    print(f"  load_taxon_hierarchy: {len(is_a)} is_a, {len(disjoint)} disjoint_from, {len(union)} union_of")


def main():
    print("test_load_constraints"); test_load_constraints()
    print("test_load_go_hierarchy"); test_load_go_hierarchy()
    print("test_load_taxon_hierarchy"); test_load_taxon_hierarchy()
    print("All data loader smoke tests passed.")


if __name__ == "__main__":
    main()
