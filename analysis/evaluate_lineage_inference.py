"""
Evaluate joint taxonomic lineage inference (no-taxons mode) across predictors
and sub-ontologies (CC, MF, BP).

Addresses Reviewer 1 Comment 3:
- Extracts both the raw CP-SAT taxon assignment (T_solver) and the GO-supported
  inferred lineage (T_forced, i.e. clades required by surviving active GO terms
  via only_in_taxon and upward taxonomy propagation).
- Maps the true NCBI Taxonomy lineage and the inferred lineage to the 7 standard
  taxonomic ranks:
    [domain (superkingdom), phylum, class, order, family, genus, species]
- Computes per (predictor, sub-ontology):
    1. Deepest correctly inferred rank (and deepest rank reachable in the GO
       taxon-constraint graph for that organism).
    2. Rank of first divergence (where an inferred clade at rank r differs from
       the true clade at rank r, vs. 'None' when the inferred lineage is an exact
       ancestral prefix of the true lineage).
    3. Cross-sub-ontology agreement across CC, MF, and BP for the same organism.
    4. Joint multi-aspect (CC+MF+BP) lineage inference accuracy when all three
       sub-ontologies share a single taxon_t vector.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict, Counter

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from taxon_consistency.adjust_ortools import (
    load_constraints,
    load_go_hierarchy,
    load_taxon_hierarchy,
    load_predictions,
    get_annotated_predictions,
    compute_demotion_flip_cost,
    solve_sat_ortools_taxon_assignment,
    get_all_ancestors_of_list,
)

RANKS_7 = ["domain", "phylum", "class", "order", "family", "genus", "species"]
DOMAIN_IDS = {
    "NCBITaxon_2": "Bacteria",
    "NCBITaxon_2157": "Archaea",
    "NCBITaxon_2759": "Eukaryota",
    "NCBITaxon_10239": "Viruses",
}


def parse_ncbitaxon_obo(obo_path: str):
    """Stream-parse ncbitaxon.obo for parent map, rank map, and name map."""
    parent_map = {}
    rank_map = {}
    name_map = {}

    cur_id = None
    cur_name = None
    cur_parent = None
    cur_rank = None

    with open(obo_path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                if cur_id:
                    if cur_name:
                        name_map[cur_id] = cur_name
                    if cur_parent:
                        parent_map[cur_id] = cur_parent
                    if cur_rank:
                        rank_map[cur_id] = cur_rank
                cur_id = cur_name = cur_parent = cur_rank = None
                continue
            if line == "[Term]":
                cur_id = cur_name = cur_parent = cur_rank = None
            elif line.startswith("id: NCBITaxon:"):
                num = line.split("id: NCBITaxon:", 1)[1].strip()
                cur_id = f"NCBITaxon_{num}"
            elif line.startswith("name: ") and cur_id:
                cur_name = line[6:].strip()
            elif line.startswith("is_a: NCBITaxon:") and cur_id:
                pnum = line.split("is_a: NCBITaxon:", 1)[1].split()[0].strip()
                cur_parent = f"NCBITaxon_{pnum}"
            elif line.startswith("property_value: has_rank NCBITaxon:") and cur_id:
                rname = line.split("property_value: has_rank NCBITaxon:", 1)[1].split()[0].strip()
                if rname in ("superkingdom", "domain", "acellular_root", "realm"):
                    cur_rank = "domain"
                elif rname in RANKS_7:
                    cur_rank = rname
                elif rname == "kingdom":
                    cur_rank = "kingdom"

        if cur_id:
            if cur_name:
                name_map[cur_id] = cur_name
            if cur_parent:
                parent_map[cur_id] = cur_parent
            if cur_rank:
                rank_map[cur_id] = cur_rank

    for did in DOMAIN_IDS:
        rank_map[did] = "domain"

    return parent_map, rank_map, name_map


def get_true_lineage(taxon_num: str, parent_map: dict):
    """Return ordered list of NCBITaxon_X from species/strain up to root."""
    cur = f"NCBITaxon_{taxon_num}"
    lineage = []
    seen = set()
    while cur and cur not in seen:
        lineage.append(cur)
        seen.add(cur)
        cur = parent_map.get(cur)
    return lineage


def lineage_by_rank(taxons_iter, rank_map: dict):
    """Map a collection of NCBITaxon_X nodes to {rank: set(taxon_ids)}."""
    by_rank = defaultdict(set)
    for t in taxons_iter:
        if t.startswith("NCBITaxon_Union_"):
            continue
        r = rank_map.get(t)
        if r in RANKS_7:
            by_rank[r].add(t)
    return by_rank


def compute_forced_taxons(
    adjusted_predictions: dict,
    in_taxon_constraints: dict,
    taxon_is_a_hierarchy: dict,
    taxon_union_of_hierarchy: dict,
    predicted_taxons: list,
):
    """
    Compute the subset of predicted_taxons that are positively supported by
    surviving GO predictions via only_in_taxon (and upward ancestor closure,
    resolving union nodes to the solver's chosen member).
    """
    active_go = set()
    for prot, tdict in adjusted_predictions.items():
        for go_id, keep in tdict.items():
            if keep:
                active_go.add(go_id)

    pred_set = set(predicted_taxons)
    seed_taxons = set()
    for go_id in active_go:
        for t in in_taxon_constraints.get(go_id, []):
            if t.startswith("NCBITaxon_Union_"):
                # Pick the union member(s) that the solver activated
                members = taxon_union_of_hierarchy.get(t, set())
                chosen = members & pred_set
                seed_taxons.update(chosen if chosen else members)
            else:
                seed_taxons.add(t)

    forced = set(seed_taxons)
    forced.update(get_all_ancestors_of_list(taxon_is_a_hierarchy, list(seed_taxons)))
    return {t for t in forced if not t.startswith("NCBITaxon_Union_")}


def evaluate_single_lineage(true_lineage: list, inferred_taxons: set, rank_map: dict, parent_map: dict):
    """
    Compare inferred_taxons against true_lineage across the 7 standard ranks.
    First expands inferred_taxons upward in NCBI Taxonomy so any inferred clade
    automatically populates its ancestor ranks (e.g. inferring Mammalia at class
    implies Chordata at phylum and Eukaryota at domain).
    """
    expanded_inferred = set()
    for t in inferred_taxons:
        if t.startswith("NCBITaxon_Union_"):
            continue
        cur = t
        seen = set()
        while cur and cur not in seen:
            expanded_inferred.add(cur)
            seen.add(cur)
            cur = parent_map.get(cur)

    true_set = set(true_lineage)
    true_ranks = {}
    for t in true_lineage:
        r = rank_map.get(t)
        if r in RANKS_7 and r not in true_ranks:
            true_ranks[r] = t

    inf_ranks = lineage_by_rank(expanded_inferred, rank_map)

    deepest_correct = "none"
    first_divergence = "none"

    for r in RANKS_7:
        t_true = true_ranks.get(r)
        t_inf_set = inf_ranks.get(r, set())
        if not t_inf_set:
            continue
        if t_true and t_inf_set == {t_true}:
            deepest_correct = r
        else:
            if first_divergence == "none":
                first_divergence = r
            break

    is_exact_prefix = expanded_inferred.issubset(true_set) and len(expanded_inferred) > 0
    return {
        "deepest_correct_rank": deepest_correct,
        "first_divergence_rank": first_divergence,
        "is_exact_prefix": is_exact_prefix,
        "true_ranks": true_ranks,
        "inferred_ranks": {k: sorted(v) for k, v in inf_ranks.items()},
        "n_inferred_nodes": len(inferred_taxons),
    }


def load_annot_and_costs_above_tau(pfile: str, tau: float):
    """Stream prediction file keeping only sc > tau to build genome-wide active GO set and flip_costs."""
    active_terms = set()
    flip_costs = defaultdict(float)
    with open(pfile, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go_id, sc_str = tok.split("|", 1)
                sc = float(sc_str)
                if sc > tau:
                    active_terms.add(go_id)
                    flip_costs[go_id] += (sc - tau)
    return {"genome": active_terms}, dict(flip_costs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--code_dir", required=True)
    ap.add_argument("--timeset_taxons_file", required=True)
    ap.add_argument("--output_json", required=True)
    ap.add_argument("--output_tsv", required=True)
    args = ap.parse_args()

    obo_path = os.path.join(args.code_dir, "data", "ncbitaxon.obo")
    constraints_file = os.path.join(args.code_dir, "data", "go_taxon_constraints_extracted_obo.tsv")
    go_hierarchy_file = os.path.join(args.code_dir, "data", "go_hierarchy.tsv")
    taxon_hierarchy_file = os.path.join(args.code_dir, "data", "taxon_hierarchy.tsv")
    ncbitaxon_hierarchy_file = os.path.join(args.code_dir, "data", "ncbitaxon_hierarchy.tsv")

    print("Parsing ncbitaxon.obo...")
    parent_map, rank_map, name_map = parse_ncbitaxon_obo(obo_path)
    print(f"Loaded {len(parent_map)} parent links, {len(rank_map)} ranked taxa.")

    in_taxon_constraints, never_in_taxon_constraints = load_constraints(constraints_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file)
    taxon_is_a, taxon_disjoint, taxon_union = load_taxon_hierarchy(
        taxon_hierarchy_file, ncbitaxon_hierarchy_file
    )

    # Compute which taxons exist in the GO taxon constraint graph
    constraint_seed = set()
    for tl in in_taxon_constraints.values():
        constraint_seed.update(tl)
    for tl in never_in_taxon_constraints.values():
        constraint_seed.update(tl)
    for u, members in taxon_union.items():
        constraint_seed.update(members)

    with open(args.timeset_taxons_file) as fh:
        timeset_taxons = sorted({line.strip() for line in fh if line.strip()})

    predictors = [
        ("seq_sim_timeset", 0.001),
        ("deepgo-se", 0.3),
        ("sprof-go", 0.001),
    ]
    subontologies = ["cc", "mf", "bp"]

    results = {}
    per_organism_records = []

    for pred_name, tau in predictors:
        results[pred_name] = {}
        # Store per-organism inferred sets across subontologies to check agreement
        org_sub_forced = defaultdict(dict)
        org_sub_solver = defaultdict(dict)
        org_joint_annot = defaultdict(dict)
        org_joint_costs = defaultdict(lambda: defaultdict(float))

        for sub in subontologies:
            pred_dir = os.path.join(
                args.data_dir, "swissprot_proteomes_folds", f"{pred_name}_results", sub, "predictions"
            )
            forced_deepest_counts = Counter()
            forced_div_counts = Counter()
            solver_deepest_counts = Counter()
            solver_div_counts = Counter()
            max_reachable_counts = Counter()
            n_exact_prefix_forced = 0
            n_exact_prefix_solver = 0
            n_evaluated = 0

            for taxon in timeset_taxons:
                matches = sorted(glob.glob(os.path.join(pred_dir, f"predictions_fold_*_taxon_{taxon}.tsv")))
                if not matches:
                    continue
                pfile = matches[0]
                annot_preds, flip_costs = load_annot_and_costs_above_tau(pfile, tau)
                # Accumulate for joint CC+MF+BP solve
                for prot, tset in annot_preds.items():
                    org_joint_annot[taxon][f"{sub}:{prot}"] = tset
                for go_id, cval in flip_costs.items():
                    org_joint_costs[taxon][go_id] += cval

                adj_preds, pred_taxons, _ = solve_sat_ortools_taxon_assignment(
                    annot_preds,
                    in_taxon_constraints,
                    never_in_taxon_constraints,
                    go_hierarchy,
                    flip_costs,
                    taxon_is_a,
                    taxon_disjoint,
                    taxon_union,
                    genome_taxon_id=None,
                )

                forced_taxons = compute_forced_taxons(
                    adj_preds, in_taxon_constraints, taxon_is_a, taxon_union, pred_taxons
                )
                solver_taxons = {t for t in pred_taxons if not t.startswith("NCBITaxon_Union_")}

                true_lin = get_true_lineage(taxon, parent_map)
                # Max reachable rank in GO constraints for this organism
                reachable_eval = evaluate_single_lineage(
                    true_lin, set(true_lin) & constraint_seed, rank_map, parent_map
                )
                max_reachable_counts[reachable_eval["deepest_correct_rank"]] += 1

                ev_forced = evaluate_single_lineage(true_lin, forced_taxons, rank_map, parent_map)
                ev_solver = evaluate_single_lineage(true_lin, solver_taxons, rank_map, parent_map)

                forced_deepest_counts[ev_forced["deepest_correct_rank"]] += 1
                forced_div_counts[ev_forced["first_divergence_rank"]] += 1
                if ev_forced["is_exact_prefix"]:
                    n_exact_prefix_forced += 1

                solver_deepest_counts[ev_solver["deepest_correct_rank"]] += 1
                solver_div_counts[ev_solver["first_divergence_rank"]] += 1
                if ev_solver["is_exact_prefix"]:
                    n_exact_prefix_solver += 1

                n_evaluated += 1
                org_sub_forced[taxon][sub] = forced_taxons
                org_sub_solver[taxon][sub] = solver_taxons

                per_organism_records.append({
                    "predictor": pred_name,
                    "subontology": sub,
                    "taxon_id": taxon,
                    "taxon_name": name_map.get(f"NCBITaxon_{taxon}", taxon),
                    "max_reachable_rank": reachable_eval["deepest_correct_rank"],
                    "forced_deepest_correct": ev_forced["deepest_correct_rank"],
                    "forced_first_divergence": ev_forced["first_divergence_rank"],
                    "forced_exact_prefix": ev_forced["is_exact_prefix"],
                    "solver_deepest_correct": ev_solver["deepest_correct_rank"],
                    "solver_first_divergence": ev_solver["first_divergence_rank"],
                    "solver_exact_prefix": ev_solver["is_exact_prefix"],
                })

            results[pred_name][sub] = {
                "n_organisms": n_evaluated,
                "max_reachable_rank_counts": dict(max_reachable_counts),
                "forced_deepest_correct_counts": dict(forced_deepest_counts),
                "forced_first_divergence_counts": dict(forced_div_counts),
                "forced_exact_prefix_count": n_exact_prefix_forced,
                "solver_deepest_correct_counts": dict(solver_deepest_counts),
                "solver_first_divergence_counts": dict(solver_div_counts),
                "solver_exact_prefix_count": n_exact_prefix_solver,
            }

        # Evaluate joint CC+MF+BP lineage inference and cross-aspect compatibility
        joint_forced_deepest = Counter()
        joint_forced_div = Counter()
        joint_solver_deepest = Counter()
        joint_solver_div = Counter()
        n_joint_exact_forced = 0
        n_joint_exact_solver = 0
        n_cross_aspect_compatible_forced = 0
        n_cross_aspect_identical_solver = 0

        for taxon in timeset_taxons:
            # Check cross-aspect compatibility of CC, MF, BP
            true_lin = get_true_lineage(taxon, parent_map)
            true_set = set(true_lin)

            # Expand each subontology's forced set in NCBI taxonomy
            sub_ranks = {}
            for sub in subontologies:
                exp = set()
                for t in org_sub_forced[taxon].get(sub, set()):
                    cur = t
                    seen = set()
                    while cur and cur not in seen:
                        exp.add(cur)
                        seen.add(cur)
                        cur = parent_map.get(cur)
                sub_ranks[sub] = lineage_by_rank(exp, rank_map)

            compatible = True
            for r in RANKS_7:
                vals = set()
                for sub in subontologies:
                    vals.update(sub_ranks[sub].get(r, set()))
                if len(vals) > 1:
                    compatible = False
                    break
            if compatible:
                n_cross_aspect_compatible_forced += 1

            # Check solver rank identity across CC, MF, BP
            solver_sub_ranks = {}
            for sub in subontologies:
                exp = set()
                for t in org_sub_solver[taxon].get(sub, set()):
                    cur = t
                    seen = set()
                    while cur and cur not in seen:
                        exp.add(cur)
                        seen.add(cur)
                        cur = parent_map.get(cur)
                solver_sub_ranks[sub] = lineage_by_rank(exp, rank_map)
            solver_compat = True
            for r in RANKS_7:
                vals = set()
                for sub in subontologies:
                    vals.update(solver_sub_ranks[sub].get(r, set()))
                if len(vals) > 1:
                    solver_compat = False
                    break
            if solver_compat:
                n_cross_aspect_identical_solver += 1

            # Joint CC+MF+BP solve
            jannot = org_joint_annot[taxon]
            jcosts = dict(org_joint_costs[taxon])
            jadj, jpred_taxons, _ = solve_sat_ortools_taxon_assignment(
                jannot,
                in_taxon_constraints,
                never_in_taxon_constraints,
                go_hierarchy,
                jcosts,
                taxon_is_a,
                taxon_disjoint,
                taxon_union,
                genome_taxon_id=None,
            )
            jforced = compute_forced_taxons(jadj, in_taxon_constraints, taxon_is_a, taxon_union, jpred_taxons)
            jsolver = {t for t in jpred_taxons if not t.startswith("NCBITaxon_Union_")}

            ev_jf = evaluate_single_lineage(true_lin, jforced, rank_map, parent_map)
            ev_js = evaluate_single_lineage(true_lin, jsolver, rank_map, parent_map)

            joint_forced_deepest[ev_jf["deepest_correct_rank"]] += 1
            joint_forced_div[ev_jf["first_divergence_rank"]] += 1
            if ev_jf["is_exact_prefix"]:
                n_joint_exact_forced += 1

            joint_solver_deepest[ev_js["deepest_correct_rank"]] += 1
            joint_solver_div[ev_js["first_divergence_rank"]] += 1
            if ev_js["is_exact_prefix"]:
                n_joint_exact_solver += 1

        results[pred_name]["joint_all_3"] = {
            "n_organisms": len(timeset_taxons),
            "cross_aspect_compatible_forced": n_cross_aspect_compatible_forced,
            "cross_aspect_compatible_solver": n_cross_aspect_identical_solver,
            "forced_deepest_correct_counts": dict(joint_forced_deepest),
            "forced_first_divergence_counts": dict(joint_forced_div),
            "forced_exact_prefix_count": n_joint_exact_forced,
            "solver_deepest_correct_counts": dict(joint_solver_deepest),
            "solver_first_divergence_counts": dict(joint_solver_div),
            "solver_exact_prefix_count": n_joint_exact_solver,
        }

    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
    with open(args.output_json, "w") as fh:
        json.dump(results, fh, indent=2)

    with open(args.output_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "predictor",
                "subontology",
                "taxon_id",
                "taxon_name",
                "max_reachable_rank",
                "forced_deepest_correct",
                "forced_first_divergence",
                "forced_exact_prefix",
                "solver_deepest_correct",
                "solver_first_divergence",
                "solver_exact_prefix",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(per_organism_records)

    print("Saved lineage inference evaluation to", args.output_json)


if __name__ == "__main__":
    main()
