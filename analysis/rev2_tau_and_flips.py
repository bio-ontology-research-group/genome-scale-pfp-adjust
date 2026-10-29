"""
Unified flip-rate and extended tau-sweep evaluation (Reviewer 1 Comments 1 & 2).

Computes for a given (predictor, subontology, tau, mode, split_scope):
1. Active predictions above tau: N_active(tau) = sum_g sum_p |{c : s(p,c) > tau}|
2. Solver grid size:             N_grid(tau)   = sum_g |P_g(tau)| * |C_g(tau)|
3. Flips above tau:              F(tau)        = sum_g sum_p |{c : s(p,c) > tau and adjusted(p,c) == False}|
4. Active-prediction flip rate:  100 * F(tau) / N_active(tau)
5. Decision-grid flip rate:      100 * F(tau) / N_grid(tau)
6. GAEF taxon consistency across all organisms in the split (evaluated at threshold tau)
7. CAFA protein-centric metrics (Fmax, Smin, AUC, AUPR) across organisms with ground-truth annotations
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "pipeline")))

from taxon_consistency.adjust_ortools import (
    load_constraints,
    load_go_hierarchy,
    load_taxon_hierarchy,
    _load_ic_map,
    solve_sat_ortools_hierarchy,
    get_all_ancestors_of_list,
)
from analysis.evaluation_scope import load_temporal_proteins, select_annotations, scope_manifest

from genome_scores_evaluator import (
    Ontology,
    load_annotations,
    evaluate_prediction_metrics,
    check_consistency_fast,
    get_ancestors_with_ontology,
    NAMESPACES,
)

FN = re.compile(r"predictions_fold_(\d+)_taxon_(\d+)\.tsv$")


def build_go_descendants(go_hierarchy: dict):
    ancestors_map = {}
    all_terms = set(go_hierarchy.keys())
    for c, ps in go_hierarchy.items():
        all_terms.update(ps)
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


def propagate_scores_dict(raw_preds: dict, go: Ontology, ont: str, annotated_proteins: set):
    namespace_set = NAMESPACES[ont]
    ancestor_cache = {}
    result = defaultdict(dict)
    for prot_id, term_scores in raw_preds.items():
        if annotated_proteins is not None and prot_id not in annotated_proteins:
            continue
        for go_id, score in term_scores.items():
            if score <= 0:
                continue
            if not go.has_term(go_id) or go.get_namespace(go_id) != namespace_set:
                continue
            if go_id not in ancestor_cache:
                ancestor_cache[go_id] = go.get_ancestors(go_id)
            if score > result[prot_id].get(go_id, 0.0):
                result[prot_id][go_id] = score
            for anc in ancestor_cache[go_id]:
                if score > result[prot_id].get(anc, 0.0):
                    result[prot_id][anc] = score
    return dict(result)


def stream_prediction_file(pfile: str, tau: float, cost_mode: str, gt_prot_ids: set):
    """
    Stream a prediction TSV file once without holding sub-threshold scores for
    all proteins in memory.
    Returns:
      annot_preds: dict[prot -> set[go_id]] for sc > tau
      flip_costs:  dict[go_id -> float]
      n_annot:     int (total predictions with sc > tau)
      raw_term_counts: dict[go_id -> int] (all raw entries per GO term)
      gt_raw_preds: dict[prot -> dict[go_id -> float]] (only for prot in gt_prot_ids)
    """
    ic_map = _load_ic_map() if cost_mode == "ic" else None
    above_tau_term_counts = defaultdict(int)
    flip_costs = defaultdict(float)
    raw_term_counts = defaultdict(int)
    gt_raw_preds = {}
    n_annot = 0
    active_prots = 0

    with open(pfile, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            prot = parts[0]
            keep_gt = gt_prot_ids is not None and prot in gt_prot_ids
            gt_dict = {} if keep_gt else None
            has_above_tau = False
            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go_id, sc_str = tok.split("|", 1)
                sc = float(sc_str)
                raw_term_counts[go_id] += 1
                if keep_gt:
                    gt_dict[go_id] = sc
                if sc > tau:
                    n_annot += 1
                    has_above_tau = True
                    above_tau_term_counts[go_id] += 1
                    if cost_mode == "uniform":
                        flip_costs[go_id] += 1.0
                    elif cost_mode == "ic":
                        flip_costs[go_id] += ic_map.get(go_id, 1.0)
                    else:
                        flip_costs[go_id] += sc - tau
            if has_above_tau:
                active_prots += 1
            if keep_gt:
                gt_raw_preds[prot] = gt_dict

    return above_tau_term_counts, active_prots, flip_costs, n_annot, raw_term_counts, gt_raw_preds


def stream_write_filtered_predictions(pfile: str, out_file: str, flipped_with_desc: set):
    """Stream pfile to out_file removing tokens whose GO term is in flipped_with_desc."""
    tmp_file = f"{out_file}.tmp.{os.getpid()}"
    with open(pfile, "r") as fin, open(tmp_file, "w") as fout:
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            kept = [parts[0]]
            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go_id = tok.split("|", 1)[0]
                if go_id not in flipped_with_desc:
                    kept.append(tok)
            fout.write("\t".join(kept) + "\n")
    os.replace(tmp_file, out_file)


def run_single_setting(
    pred_name: str,
    sub: str,
    tau: float,
    provide_taxon_id: bool,
    cost_mode: str,
    target_taxons: set,
    scope_label: str,
    data_dir: str,
    code_dir: str,
    go: Ontology,
    in_taxon: dict,
    never_in_taxon: dict,
    go_hierarchy: dict,
    go_descendants: dict,
    taxon_is_a: dict,
    taxon_disjoint: dict,
    taxon_union: dict,
    eval_cafa: bool = True,
    save_dir: str = None,
    temporal_proteins: dict = None,
):
    os.environ["PFP_COST_MODE"] = cost_mode
    pred_dir = os.path.join(
        data_dir, "swissprot_proteomes_folds", f"{pred_name}_results", sub, "predictions"
    )
    ann_dir = os.path.join(data_dir, "swissprot_proteomes_folds", "annotations-go-basic")
    taxa_constraints_file = os.path.join(code_dir, "data", "constraints", "taxon_constraints.tsv")

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    total_active_preds = 0
    total_grid_decisions = 0
    total_flips_above_tau = 0
    total_raw_removed = 0
    n_organisms = 0

    gaef_consistent_before = []
    gaef_consistent_after = []

    cafa_fmax = []
    cafa_smin = []
    cafa_auc = []
    cafa_aupr = []
    cafa_wfmax = []

    pooled_annotations = {}
    pooled_predictions = {}

    for pfile in sorted(glob.glob(os.path.join(pred_dir, "predictions_fold_*_taxon_*.tsv"))):
        m = FN.search(os.path.basename(pfile))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        if target_taxons is not None and taxon not in target_taxons:
            continue

        n_organisms += 1
        org_annots = None
        gt_prot_ids = None
        if eval_cafa:
            ann_file = os.path.join(ann_dir, f"annots_taxon_{taxon}.tsv")
            if os.path.exists(ann_file):
                org_annots = select_annotations(load_annotations(ann_file), taxon, temporal_proteins)
                if org_annots:
                    gt_prot_ids = set(org_annots.keys())

        above_tau_term_counts, active_prots, flip_costs, n_annot, raw_term_counts, gt_raw_preds = stream_prediction_file(
            pfile, tau, cost_mode, gt_prot_ids
        )

        # Grid size |P_g(tau)| * |C_g(tau)|
        active_terms = set(above_tau_term_counts.keys())
        grid_size = active_prots * len(active_terms)

        adj_preds, _ = solve_sat_ortools_hierarchy(
            {"genome": active_terms},
            in_taxon,
            never_in_taxon,
            go_hierarchy,
            flip_costs,
            taxon_is_a,
            taxon_disjoint,
            taxon_union,
            genome_taxon_id=taxon if provide_taxon_id else None,
        )

        # Determine flipped GO terms (and their descendants)
        flipped_terms = set()
        for p, tdict in adj_preds.items():
            for g, keep in tdict.items():
                if not keep:
                    flipped_terms.add(g)
        del adj_preds
        n_flips = sum(above_tau_term_counts.get(g, 0) for g in flipped_terms)
        del above_tau_term_counts
        flipped_with_desc = set(flipped_terms)
        for g in flipped_terms:
            flipped_with_desc.update(go_descendants.get(g, ()))

        raw_removed = sum(raw_term_counts.get(g, 0) for g in flipped_with_desc)
        del raw_term_counts

        if save_dir:
            out_file = os.path.join(save_dir, f"optimized_fold_{fold}_taxon_{taxon}.tsv")
            if not os.path.exists(out_file):
                stream_write_filtered_predictions(pfile, out_file, flipped_with_desc)

        total_active_preds += n_annot
        total_grid_decisions += grid_size
        total_flips_above_tau += n_flips
        total_raw_removed += raw_removed

        # GAEF taxonomic consistency at threshold tau (depends only on genome-wide union of active GO terms)
        before_anc, _ = get_ancestors_with_ontology({"genome": active_terms}, go)
        after_anc, _ = get_ancestors_with_ontology({"genome": active_terms - flipped_with_desc}, go)
        gaef_consistent_before.append(bool(check_consistency_fast(before_anc, taxa_constraints_file)))
        gaef_consistent_after.append(bool(check_consistency_fast(after_anc, taxa_constraints_file)))

        if eval_cafa and org_annots:
            gt_final_preds = {
                p: {g: sc for g, sc in gdict.items() if g not in flipped_with_desc}
                for p, gdict in gt_raw_preds.items()
            }
            prop_preds = propagate_scores_dict(gt_final_preds, go, sub, gt_prot_ids)
            m_org = evaluate_prediction_metrics(go, sub, org_annots, prop_preds)
            if m_org is not None:
                cafa_fmax.append(m_org["fmax"])
                if np.isfinite(m_org["smin"]):
                    cafa_smin.append(m_org["smin"])
                cafa_auc.append(m_org["avg_auc"])
                cafa_aupr.append(m_org["aupr"])
                cafa_wfmax.append(m_org["wfmax"])
                for prot_id, gt_terms in org_annots.items():
                    pooled_annotations[f"{taxon}:{prot_id}"] = gt_terms
                    if prot_id in prop_preds:
                        pooled_predictions[f"{taxon}:{prot_id}"] = prop_preds[prot_id]

    pooled_metrics = {}
    if eval_cafa and pooled_annotations:
        pm = evaluate_prediction_metrics(go, sub, pooled_annotations, pooled_predictions)
        if pm is not None:
            pooled_metrics = {
                "fmax": round(float(pm["fmax"]), 4),
                "smin": round(float(pm["smin"]), 2),
                "avg_auc": round(float(pm["avg_auc"]), 4),
                "aupr": round(float(pm["aupr"]), 4),
                "wfmax": round(float(pm["wfmax"]), 4),
                "tmax": round(float(pm["tmax"]), 2),
            }

    active_rate = (100.0 * total_flips_above_tau / total_active_preds) if total_active_preds else 0.0
    grid_rate = (100.0 * total_flips_above_tau / total_grid_decisions) if total_grid_decisions else 0.0

    return {
        "predictor": pred_name,
        "subontology": sub,
        "scope": scope_label,
        "mode": "with_taxons" if provide_taxon_id else "no_taxons",
        "cost_mode": cost_mode,
        "tau": tau,
        "n_organisms": n_organisms,
        "n_cafa_organisms": len(cafa_fmax),
        "n_annotation_proteins": len(pooled_annotations),
        "evaluation_scope": "entry_date_mask" if temporal_proteins is not None else "organism_holdout",
        "preds_above_tau": total_active_preds,
        "grid_decisions": total_grid_decisions,
        "flips_above_tau": total_flips_above_tau,
        "active_flip_rate_pct": round(active_rate, 4),
        "grid_flip_rate_pct": round(grid_rate, 4),
        "raw_entries_removed": total_raw_removed,
        "gaef_taxon_consistency_before": round(float(np.mean(gaef_consistent_before)), 4) if gaef_consistent_before else None,
        "gaef_taxon_consistency_after": round(float(np.mean(gaef_consistent_after)), 4) if gaef_consistent_after else None,
        "macro_cafa": {
            "fmax": round(float(np.mean(cafa_fmax)), 4) if cafa_fmax else None,
            "smin": round(float(np.mean(cafa_smin)), 2) if cafa_smin else None,
            "avg_auc": round(float(np.mean(cafa_auc)), 4) if cafa_auc else None,
            "aupr": round(float(np.mean(cafa_aupr)), 4) if cafa_aupr else None,
        },
        "pooled_cafa": pooled_metrics,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--code_dir", required=True)
    ap.add_argument("--timeset_taxons_file", required=True)
    ap.add_argument("--task", required=True, choices=["table1_flips", "tau_sweep", "cost_ablation"])
    ap.add_argument("--predictor", default="deepgo-se")
    ap.add_argument("--subontology", default="cc")
    ap.add_argument("--output_json", required=True)
    ap.add_argument("--save_adjusted_root", default=None)
    ap.add_argument("--temporal-proteins", help="Explicit entry-date protein TSV; required for timeset accuracy")
    args = ap.parse_args()
    if args.predictor in ("seq_sim_timeset", "deepgo-se", "sprof-go") and not args.temporal_proteins:
        ap.error("--temporal-proteins is required: organism membership is not a protein test mask")
    temporal_proteins = load_temporal_proteins(args.temporal_proteins) if args.temporal_proteins else None

    with open(args.timeset_taxons_file) as fh:
        timeset_taxons = {line.strip() for line in fh if line.strip()}

    constraints_file = os.path.join(args.code_dir, "data", "go_taxon_constraints_extracted_obo.tsv")
    go_hierarchy_file = os.path.join(args.code_dir, "data", "go_hierarchy.tsv")
    taxon_hierarchy_file = os.path.join(args.code_dir, "data", "taxon_hierarchy.tsv")
    ncbitaxon_hierarchy_file = os.path.join(args.code_dir, "data", "ncbitaxon_hierarchy.tsv")
    go_obo_file = os.path.join(args.code_dir, "data", "go-basic.obo")

    in_taxon, never_in_taxon = load_constraints(constraints_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file)
    go = Ontology(go_obo_file)
    for t_id, t_info in go.ont.items():
        for rel in ("is_a", "part_of"):
            for p_id in t_info.get(rel, ()):
                go_hierarchy[t_id].add(p_id)
    _, go_descendants = build_go_descendants(go_hierarchy)
    taxon_is_a, taxon_disjoint, taxon_union = load_taxon_hierarchy(taxon_hierarchy_file, ncbitaxon_hierarchy_file)

    default_taus = {
        "seq_sim_timeset": 0.001,
        "deepgo-se": 0.3,
        "sprof-go": 0.001,
        "seq_sim_heldout": 0.001,
        "mlp_heldout": 0.001,
    }

    records = []
    if args.task == "table1_flips":
        tau = default_taus[args.predictor]
        is_timeset_pred = args.predictor in ("seq_sim_timeset", "deepgo-se", "sprof-go")
        scopes = [("timeset_50", timeset_taxons)] if is_timeset_pred else [("heldout_112", None)]
        if args.predictor == "deepgo-se" and args.subontology == "cc":
            scopes.append(("all_564", None))

        for scope_label, tset in scopes:
            for provide_tax in [True, False]:
                # Only run expensive no-taxons on the test split (timeset_50 or heldout_112)
                if scope_label == "all_564" and not provide_tax:
                    continue
                sdir = None
                if args.save_adjusted_root and scope_label == "timeset_50" and args.subontology == "cc":
                    mname = "optimized" if provide_tax else "optimized_without_taxons"
                    sdir = os.path.join(args.save_adjusted_root, args.predictor, args.subontology, mname)
                rec = run_single_setting(
                    args.predictor,
                    args.subontology,
                    tau,
                    provide_tax,
                    "margin",
                    tset,
                    scope_label,
                    args.data_dir,
                    args.code_dir,
                    go,
                    in_taxon,
                    never_in_taxon,
                    go_hierarchy,
                    go_descendants,
                    taxon_is_a,
                    taxon_disjoint,
                    taxon_union,
                    eval_cafa=(scope_label in ("timeset_50", "heldout_112") and args.subontology != "bp"),
                    save_dir=sdir,
                    temporal_proteins=temporal_proteins,
                )
                records.append(rec)
                print(json.dumps(rec))

    elif args.task == "tau_sweep":
        is_timeset_pred = args.predictor in ("seq_sim_timeset", "deepgo-se", "sprof-go")
        tset = timeset_taxons if is_timeset_pred else None
        scope_label = "timeset_50" if is_timeset_pred else "heldout_112"
        for tau in [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
            rec = run_single_setting(
                args.predictor,
                args.subontology,
                tau,
                True,
                "margin",
                tset,
                scope_label,
                args.data_dir,
                args.code_dir,
                go,
                in_taxon,
                never_in_taxon,
                go_hierarchy,
                go_descendants,
                taxon_is_a,
                taxon_disjoint,
                taxon_union,
                eval_cafa=True,
                save_dir=None,
                temporal_proteins=temporal_proteins,
            )
            records.append(rec)
            print(json.dumps(rec))

    elif args.task == "cost_ablation":
        os.environ.setdefault(
            "PFP_IC_FILE",
            "/ibex/scratch/projects/c2014/hohndor/pfp-adjust-rh/out/cost_ablation/go_ic.tsv",
        )
        tau = default_taus[args.predictor]
        for scope_label, tset in [("timeset_50", timeset_taxons), ("all_564", None)]:
            for cmode in ["margin", "uniform", "ic"]:
                rec = run_single_setting(
                    args.predictor,
                    args.subontology,
                    tau,
                    True,
                    cmode,
                    tset,
                    scope_label,
                    args.data_dir,
                    args.code_dir,
                    go,
                    in_taxon,
                    never_in_taxon,
                    go_hierarchy,
                    go_descendants,
                    taxon_is_a,
                    taxon_disjoint,
                    taxon_union,
                    eval_cafa=(scope_label == "timeset_50"),
                    save_dir=None,
                temporal_proteins=temporal_proteins,
                )
                records.append(rec)
                print(json.dumps(rec))

    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)
    with open(args.output_json, "w") as fh:
        json.dump(records, fh, indent=2)


if __name__ == '__main__':
    raise SystemExit('This historical evaluator is retired. Use the documented reproducibility/revision workflow.')
