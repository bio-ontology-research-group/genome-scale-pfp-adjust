"""
HISTORICAL EXPLORATORY IMPLEMENTATION; NOT ACCEPTED REVISION EVIDENCE.
The current workflow is analysis/revision_audit_run.py.

Historical local corpus-absence filter inspired by FunTaxIS-lite.

Implements a local corpus-absence approximation, not the published FunTaxIS-lite pipeline:
1. Partitions NCBI Taxonomy into the 172 reference clades from
   data/funtaxis_lite/taxonConstraintsDef.txt ('reliable' vs 'unreliable').
2. Computes cumulative GO term annotation frequencies across the training
   proteomes (excluding the 50 timeset test organisms to prevent test leakage)
   with true-path upward propagation.
3. For every GO term with global frequency >= cutoff (default 500) that has 0
   annotations in a 'reliable' reference clade T, generates an automatic
   never_in_taxon(c, T) rule and propagates it downward to all GO descendants of c.
4. Merges the automatic empirical rules with the GO Consortium's official
   only_in_taxon and never_in_taxon constraints (propagated downward to GO
   descendants) and data/funtaxis_lite/manualConstraints.txt.
5. Filters prediction files for each test organism by removing all GO terms in
   that organism's merged never_in set.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from taxon_consistency.adjust_ortools import (
    load_constraints,
    load_go_hierarchy,
    load_taxon_hierarchy,
    get_all_ancestors_of_list,
)

FN = re.compile(r"predictions_fold_(\d+)_taxon_(\d+)\.tsv$")


def build_go_ancestors_and_descendants(go_hierarchy: dict):
    """Build full transitive ancestors and descendants maps for GO terms."""
    children_map = defaultdict(set)
    all_terms = set(go_hierarchy.keys())
    for child, parents in go_hierarchy.items():
        for p in parents:
            children_map[p].add(child)
            all_terms.add(p)

    ancestors_map = {}
    for t in all_terms:
        anc = set()
        stack = list(go_hierarchy.get(t, ()))
        while stack:
            cur = stack.pop()
            if cur not in anc:
                anc.add(cur)
                stack.extend(go_hierarchy.get(cur, ()))
        ancestors_map[t] = anc

    descendants_map = defaultdict(set)
    for t, ancs in ancestors_map.items():
        for a in ancs:
            descendants_map[a].add(t)

    return ancestors_map, descendants_map


def load_funtaxis_partitions(partition_file: str):
    """Load taxonConstraintsDef.txt -> {NCBITaxon_X: ('reliable'|'unreliable', name, rank)}."""
    partitions = {}
    with open(partition_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0].isdigit():
                tid = f"NCBITaxon_{parts[0]}"
                partitions[tid] = {
                    "reliability": parts[2].strip(),
                    "name": parts[1].strip(),
                    "rank": parts[5].strip() if len(parts) > 5 else "",
                }
    return partitions


def load_manual_constraints(manual_file: str):
    """Load manualConstraints.txt -> (manual_in, manual_never) mapping taxon -> set(go_terms)."""
    manual_in = defaultdict(set)
    manual_never = defaultdict(set)
    with open(manual_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            go_id = parts[0].replace("_", ":")
            tid = f"NCBITaxon_{parts[3].strip()}"
            ctype = parts[5].strip().lower()
            if ctype == "never in taxon":
                manual_never[tid].add(go_id)
            elif ctype == "in taxon":
                manual_in[tid].add(go_id)
    return manual_in, manual_never


def compute_funtaxis_forbidden_terms(
    target_taxons: list,
    annotations_dir: str,
    timeset_taxons: set,
    go_ancestors: dict,
    go_descendants: dict,
    in_taxon_constraints: dict,
    never_in_taxon_constraints: dict,
    taxon_is_a: dict,
    taxon_union: dict,
    partitions: dict,
    manual_in: dict,
    manual_never: dict,
    cutoff: int = 500,
):
    """
    Compute the FunTaxIS-lite forbidden GO term set for each taxon in target_taxons.
    Returns:
      forbidden_by_taxon: dict[taxon_str -> set[go_id]]
      auto_only_by_taxon: dict[taxon_str -> set[go_id]] (terms added solely by empirical rules)
      go_consortium_by_taxon: dict[taxon_str -> set[go_id]] (terms from official GO axioms)
    """
    # Step 1: Compute cumulative GO term frequencies across training organisms
    global_go_freq = defaultdict(int)
    clade_go_freq = defaultdict(lambda: defaultdict(int))

    part_nodes = set(partitions.keys())

    for ann_path in sorted(glob.glob(os.path.join(annotations_dir, "annots_taxon_*.tsv"))):
        base = os.path.basename(ann_path)
        tnum = base.replace("annots_taxon_", "").replace(".tsv", "")
        if tnum in timeset_taxons:
            # Exclude test split organisms to prevent data leakage
            continue
        tnode = f"NCBITaxon_{tnum}"
        t_ancs = get_all_ancestors_of_list(taxon_is_a, [tnode]) | {tnode}
        matching_partitions = t_ancs & part_nodes

        with open(ann_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                raw_terms = set(parts[1:])
                prop_terms = set(raw_terms)
                for g in raw_terms:
                    prop_terms.update(go_ancestors.get(g, ()))
                for g in prop_terms:
                    global_go_freq[g] += 1
                    for pnode in matching_partitions:
                        clade_go_freq[pnode][g] += 1

    eligible_go = {g for g, cnt in global_go_freq.items() if cnt >= cutoff}
    print(f"[FunTaxIS-lite] Training eligible GO terms (freq >= {cutoff}): {len(eligible_go)}")

    # Step 2: For each 'reliable' partition clade, find eligible GO terms with 0 count
    auto_never_per_partition = defaultdict(set)
    for pnode, pinfo in partitions.items():
        if pinfo["reliability"] != "reliable":
            continue
        pcounts = clade_go_freq[pnode]
        for g in eligible_go:
            if pcounts.get(g, 0) == 0:
                auto_never_per_partition[pnode].add(g)
                auto_never_per_partition[pnode].update(go_descendants.get(g, ()))

    # Expand union definitions for official GO constraints
    expanded_in = {}
    for g, tlist in in_taxon_constraints.items():
        expanded = set()
        for t in tlist:
            if t.startswith("NCBITaxon_Union_"):
                expanded.update(taxon_union.get(t, ()))
            else:
                expanded.add(t)
        expanded_in[g] = expanded

    expanded_never = {}
    for g, tlist in never_in_taxon_constraints.items():
        expanded = set()
        for t in tlist:
            if t.startswith("NCBITaxon_Union_"):
                expanded.update(taxon_union.get(t, ()))
            else:
                expanded.add(t)
        expanded_never[g] = expanded

    forbidden_by_taxon = {}
    auto_only_by_taxon = {}
    go_consortium_by_taxon = {}

    for tnum in target_taxons:
        tnode = f"NCBITaxon_{tnum}"
        t_ancs = get_all_ancestors_of_list(taxon_is_a, [tnode]) | {tnode}

        # Manual 'in taxon' overrules on this lineage (propagated to GO descendants)
        man_in_terms = set()
        for anc in t_ancs:
            for g in manual_in.get(anc, ()):
                man_in_terms.add(g)
                man_in_terms.update(go_descendants.get(g, ()))

        # Official GO Consortium allowed (only_in satisfied) terms
        go_only_in_satisfied = set()
        go_consortium_never = set()
        for g, req_taxa in expanded_in.items():
            if req_taxa & t_ancs:
                go_only_in_satisfied.add(g)
            else:
                go_consortium_never.add(g)
                go_consortium_never.update(go_descendants.get(g, ()))

        for g, forb_taxa in expanded_never.items():
            if forb_taxa & t_ancs:
                go_consortium_never.add(g)
                go_consortium_never.update(go_descendants.get(g, ()))

        go_consortium_never -= man_in_terms

        # Closest partition ancestor in taxonConstraintsDef.txt
        # Walk up t_ancs: apply automatic never_in from reliable partition ancestors
        auto_never = set()
        for anc in t_ancs & part_nodes:
            auto_never.update(auto_never_per_partition.get(anc, ()))
        auto_never -= (man_in_terms | go_only_in_satisfied)

        # Manual 'never in taxon' on this lineage
        man_never_terms = set()
        for anc in t_ancs:
            for g in manual_never.get(anc, ()):
                man_never_terms.add(g)
                man_never_terms.update(go_descendants.get(g, ()))

        merged_forbidden = go_consortium_never | auto_never | man_never_terms
        forbidden_by_taxon[tnum] = merged_forbidden
        auto_only_by_taxon[tnum] = merged_forbidden - go_consortium_never
        go_consortium_by_taxon[tnum] = go_consortium_never

    return forbidden_by_taxon, auto_only_by_taxon, go_consortium_by_taxon


def filter_predictions_dir(
    predictions_dir: str,
    output_dir: str,
    target_taxons: set,
    forbidden_by_taxon: dict,
    tau: float,
):
    """Filter prediction files for target_taxons and return flip statistics at threshold tau."""
    os.makedirs(output_dir, exist_ok=True)
    preds_above_tau = 0
    flips_above_tau = 0
    n_files = 0

    for pfile in sorted(glob.glob(os.path.join(predictions_dir, "predictions_fold_*_taxon_*.tsv"))):
        m = FN.search(os.path.basename(pfile))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        if target_taxons and taxon not in target_taxons:
            continue
        forb = forbidden_by_taxon.get(taxon, set())
        out_file = os.path.join(output_dir, f"optimized_fold_{fold}_taxon_{taxon}.tsv")
        n_files += 1

        with open(pfile) as fin, open(out_file, "w") as fout:
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0]:
                    continue
                prot = parts[0]
                kept = [prot]
                for tok in parts[1:]:
                    if "|" not in tok:
                        continue
                    go_id, sc_str = tok.split("|", 1)
                    sc = float(sc_str)
                    if sc > tau:
                        preds_above_tau += 1
                    if go_id in forb:
                        if sc > tau:
                            flips_above_tau += 1
                    else:
                        kept.append(tok)
                fout.write("\t".join(kept) + "\n")

    rate = (100.0 * flips_above_tau / preds_above_tau) if preds_above_tau else 0.0
    return {
        "n_organisms": n_files,
        "preds_above_tau": preds_above_tau,
        "flips_above_tau": flips_above_tau,
        "flip_rate_pct": round(rate, 4),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--code_dir", required=True)
    ap.add_argument("--timeset_taxons_file", required=True)
    ap.add_argument("--out_root", required=True)
    ap.add_argument("--summary_json", required=True)
    ap.add_argument("--cutoff", type=int, default=500)
    args = ap.parse_args()

    with open(args.timeset_taxons_file) as fh:
        timeset_taxons = sorted({line.strip() for line in fh if line.strip()})
    timeset_set = set(timeset_taxons)

    constraints_file = os.path.join(args.code_dir, "data", "go_taxon_constraints_extracted_obo.tsv")
    go_hierarchy_file = os.path.join(args.code_dir, "data", "go_hierarchy.tsv")
    taxon_hierarchy_file = os.path.join(args.code_dir, "data", "taxon_hierarchy.tsv")
    ncbitaxon_hierarchy_file = os.path.join(args.code_dir, "data", "ncbitaxon_hierarchy.tsv")
    partition_file = os.path.join(args.code_dir, "data", "funtaxis_lite", "taxonConstraintsDef.txt")
    manual_file = os.path.join(args.code_dir, "data", "funtaxis_lite", "manualConstraints.txt")
    annotations_dir = os.path.join(args.data_dir, "swissprot_proteomes_folds", "annotations-go-basic")

    in_taxon, never_in_taxon = load_constraints(constraints_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file)
    taxon_is_a, _, taxon_union = load_taxon_hierarchy(taxon_hierarchy_file, ncbitaxon_hierarchy_file)
    go_ancestors, go_descendants = build_go_ancestors_and_descendants(go_hierarchy)
    partitions = load_funtaxis_partitions(partition_file)
    manual_in, manual_never = load_manual_constraints(manual_file)

    forbidden_by_taxon, auto_only_by_taxon, go_cons_by_taxon = compute_funtaxis_forbidden_terms(
        timeset_taxons,
        annotations_dir,
        timeset_set,
        go_ancestors,
        go_descendants,
        in_taxon,
        never_in_taxon,
        taxon_is_a,
        taxon_union,
        partitions,
        manual_in,
        manual_never,
        cutoff=args.cutoff,
    )

    avg_total_forb = sum(len(v) for v in forbidden_by_taxon.values()) / len(timeset_taxons)
    avg_auto_only = sum(len(v) for v in auto_only_by_taxon.values()) / len(timeset_taxons)
    avg_go_cons = sum(len(v) for v in go_cons_by_taxon.values()) / len(timeset_taxons)
    print(
        f"[FunTaxIS-lite] Mean forbidden GO terms per timeset organism: "
        f"total={avg_total_forb:.1f} (GO-Consortium={avg_go_cons:.1f}, empirical-auto-extra={avg_auto_only:.1f})"
    )

    predictors = [
        ("seq_sim_timeset", 0.001),
        ("deepgo-se", 0.3),
        ("sprof-go", 0.001),
    ]
    summary = {
        "cutoff": args.cutoff,
        "mean_forbidden_terms_total": round(avg_total_forb, 1),
        "mean_forbidden_terms_go_consortium": round(avg_go_cons, 1),
        "mean_forbidden_terms_empirical_extra": round(avg_auto_only, 1),
        "runs": {},
    }

    for pred_name, tau in predictors:
        summary["runs"][pred_name] = {}
        for sub in ["cc", "mf", "bp"]:
            pdir = os.path.join(
                args.data_dir, "swissprot_proteomes_folds", f"{pred_name}_results", sub, "predictions"
            )
            odir = os.path.join(args.out_root, pred_name, sub, "funtaxis_lite")
            stats = filter_predictions_dir(pdir, odir, timeset_set, forbidden_by_taxon, tau)
            summary["runs"][pred_name][sub] = stats
            print(f"[FunTaxIS-lite] {pred_name} {sub}: {stats}")

    os.makedirs(os.path.dirname(args.summary_json), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == '__main__':
    raise SystemExit('This historical evaluator is retired. Use analysis/revision_audit_run.py and the documented revision workflow.')
