"""
Comprehensive evaluation of Stage 2 modes (standard vs demotion-only) and
external baselines (FunTaxIS-lite, Naive filter, Soft Regularizer / SBR)
across the 50 timeset organisms (Reviewer 1 Comments 5, 6, 7).

Reports:
- GAEF taxon_consistency and complex_coherence over ALL 50 timeset organisms
  AND over the 41 timeset organisms with CC ground truth (46 MF, 46 BP).
- CAFA metrics (Fmax, Smin, AUC, AUPR) both micro-pooled (matching Table 2)
  and macro-averaged across organisms.
- Stage 2 Complex Portal validation for both standard and demotion-only modes.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
import time
from collections import defaultdict
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "pipeline")))

from taxon_consistency.adjust_ortools import (
    load_constraints,
    load_go_hierarchy as load_full_go_hierarchy,
    load_taxon_hierarchy,
)
from complex_coherence.adjust_ortools import (
    main as complex_coherence_adjust,
    load_go_hierarchy as load_cc_go_hierarchy,
    load_homodimer_terms,
    get_heteromeric_complexes,
)
from genome_scores_evaluator import (
    Ontology,
    load_annotations,
    evaluate_prediction_metrics,
    check_consistency_fast,
    get_ancestors_with_ontology,
    NAMESPACES,
)
from gaef_patches import coherence
from gaef_patches.complex_classifier import parse_homodimer_terms_file

FN_PRED = re.compile(r"_fold_(\d+)_taxon_(\d+)\.tsv$")


def transitive(seed, adj):
    out, stack = set(seed), list(seed)
    while stack:
        for nxt in adj.get(stack.pop(), ()):
            if nxt not in out:
                out.add(nxt)
                stack.append(nxt)
    return out


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


_COMPLEX_ANC_CACHE = {}


def _get_complex_ancestors(go_term: str, go_is_a: Ontology, complex_child_terms: set):
    res = _COMPLEX_ANC_CACHE.get(go_term)
    if res is None:
        anc = go_is_a.get_ancestors(go_term)
        res = frozenset(anc & complex_child_terms) if "GO:0032991" in anc else frozenset()
        _COMPLEX_ANC_CACHE[go_term] = res
    return res


def fast_gaef_eval_from_unions(
    union_terms: set,
    complex_prot_counts: dict,
    sub: str,
    go: Ontology,
    taxa_constraints_file: str,
    homodimer_terms_gaef: set,
):
    """
    Compute GAEF satisfiable (taxon consistency) and complex_coherence (for CC)
    directly from the genome-wide union of predicted terms and per-complex protein counts.
    """
    anc, _ = get_ancestors_with_ontology({"genome": union_terms}, go)
    is_sat = bool(check_consistency_fast(anc, taxa_constraints_file))
    c_coh = None
    if sub == "cc":
        coh_cnt = 0
        incoh_cnt = 0
        for c, cnt in complex_prot_counts.items():
            if (homodimer_terms_gaef and c in homodimer_terms_gaef) or cnt > 1:
                coh_cnt += 1
            else:
                incoh_cnt += 1
        c_coh = (coh_cnt / (coh_cnt + incoh_cnt)) * 100.0 if (coh_cnt + incoh_cnt) > 0 else 0.0
    return is_sat, c_coh


def stream_file_for_eval(
    fpath: str,
    gaef_tau: float,
    gt_prot_ids: set,
    complex_anc_map: dict,
    forbidden_terms: set = None,
):
    """
    Stream a prediction TSV file once and return:
      union_tau: set[go_id] with sc > gaef_tau across genome
      comp_counts_tau: dict[complex_id -> int] of proteins with sc > gaef_tau
      union_03: set[go_id] with sc > 0.3 across genome
      comp_counts_03: dict[complex_id -> int] of proteins with sc > 0.3
      gt_raw_preds: dict[prot -> dict[go_id -> float]] for prot in gt_prot_ids
    """
    union_tau = set()
    comp_counts_tau = defaultdict(int)
    union_03 = set()
    comp_counts_03 = defaultdict(int)
    gt_raw_preds = {}
    has_forb = bool(forbidden_terms)
    has_comp = bool(complex_anc_map)

    with open(fpath, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            prot = parts[0]
            keep_gt = gt_prot_ids is not None and prot in gt_prot_ids
            gt_dict = {} if keep_gt else None
            prot_c_tau = set() if has_comp else None
            prot_c_03 = set() if has_comp else None

            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go_id, sc_str = tok.split("|", 1)
                if has_forb and go_id in forbidden_terms:
                    continue
                ca = complex_anc_map.get(go_id) if has_comp else None
                if not keep_gt and ca is None and (go_id in union_tau) and (go_id in union_03):
                    continue
                sc = float(sc_str)
                if sc > gaef_tau:
                    union_tau.add(go_id)
                    if ca:
                        prot_c_tau.update(ca)
                if sc > 0.3:
                    union_03.add(go_id)
                    if ca:
                        prot_c_03.update(ca)
                if keep_gt:
                    gt_dict[go_id] = sc

            if has_comp:
                for c in prot_c_tau:
                    comp_counts_tau[c] += 1
                for c in prot_c_03:
                    comp_counts_03[c] += 1
            if keep_gt:
                gt_raw_preds[prot] = gt_dict

    return union_tau, comp_counts_tau, union_03, comp_counts_03, gt_raw_preds


_WORKER_CTX = {}


def _init_worker(ctx):
    global _WORKER_CTX
    _WORKER_CTX = ctx


def _eval_single_organism_worker(task):
    taxon, fold, fpath, gaef_tau, sub, forb_mode = task
    ctx = _WORKER_CTX
    go = ctx["go"]
    taxa_constraints_file = ctx["taxa_constraints_file"]
    homodimer_terms_gaef = ctx["homodimer_terms_gaef"]
    complex_anc_map = ctx["complex_anc_map"] if sub == "cc" else None
    ann_dir = ctx["ann_dir"]
    if forb_mode == "noprop":
        forb = ctx["naive_forb_noprop"].get(taxon)
    elif forb_mode == "prop":
        forb = ctx["naive_forb_prop"].get(taxon)
    else:
        forb = None

    ann_file = os.path.join(ann_dir, f"annots_taxon_{taxon}.tsv")
    org_annots = load_annotations(ann_file) if os.path.exists(ann_file) else None
    gt_prot_ids = set(org_annots.keys()) if org_annots else None

    union_tau, comp_counts_tau, union_03, comp_counts_03, gt_raw_preds = stream_file_for_eval(
        fpath, gaef_tau, gt_prot_ids, complex_anc_map, forb
    )
    is_sat, c_coh = fast_gaef_eval_from_unions(
        union_tau, comp_counts_tau, sub, go, taxa_constraints_file, homodimer_terms_gaef
    )
    if abs(gaef_tau - 0.3) > 1e-6:
        is_sat03, c_coh03 = fast_gaef_eval_from_unions(
            union_03, comp_counts_03, sub, go, taxa_constraints_file, homodimer_terms_gaef
        )
    else:
        is_sat03, c_coh03 = is_sat, c_coh

    has_gt = False
    prop_preds = {}
    if org_annots:
        go_set = set(go.get_namespace_terms(NAMESPACES[sub]))
        go_set.discard({"cc": "GO:0005575", "mf": "GO:0003674", "bp": "GO:0008150"}[sub])
        has_sub_gt = any((terms & go_set) for terms in org_annots.values())
        if has_sub_gt:
            has_gt = True
            prop_preds = propagate_scores_dict(gt_raw_preds, go, sub, gt_prot_ids)

    return (taxon, is_sat, c_coh, is_sat03, c_coh03, has_gt, org_annots if has_gt else None, prop_preds if has_gt else None)


def evaluate_dir_on_timeset(
    dir_path: str,
    sub: str,
    gaef_tau: float,
    timeset_taxons: set,
    ann_dir: str,
    go: Ontology,
    go_is_a: Ontology,
    taxa_constraints_file: str,
    complex_child_terms: set,
    homodimer_terms_gaef: set,
    forbidden_by_taxon: dict = None,
    pool=None,
):
    tc_all50 = []
    cc_all50 = []
    tc_gt = []
    cc_gt = []

    tc_all50_03 = []
    cc_all50_03 = []
    tc_gt_03 = []
    cc_gt_03 = []

    pooled_annotations = {}
    pooled_predictions = {}

    files_by_taxon = {}
    for fpath in sorted(glob.glob(os.path.join(dir_path, "*_fold_*_taxon_*.tsv"))):
        m = FN_PRED.search(os.path.basename(fpath))
        if m:
            fold, taxon = m.group(1), m.group(2)
            if taxon in timeset_taxons:
                files_by_taxon[taxon] = (fold, fpath)

    tasks = []
    for taxon in sorted(files_by_taxon.keys()):
        fold, fpath = files_by_taxon[taxon]
        tasks.append((taxon, fold, fpath, gaef_tau, sub, forbidden_by_taxon))

    results = list(pool.map(_eval_single_organism_worker, tasks)) if pool is not None else [
        _eval_single_organism_worker(t) for t in tasks
    ]

    for taxon, is_sat, c_coh, is_sat03, c_coh03, has_gt, org_annots, prop_preds in results:
        tc_all50.append(is_sat)
        if c_coh is not None:
            cc_all50.append(c_coh)
        tc_all50_03.append(is_sat03)
        if c_coh03 is not None:
            cc_all50_03.append(c_coh03)

        if has_gt:
            tc_gt.append(is_sat)
            if c_coh is not None:
                cc_gt.append(c_coh)
            tc_gt_03.append(is_sat03)
            if c_coh03 is not None:
                cc_gt_03.append(c_coh03)
            for prot_id, gt_terms in org_annots.items():
                pooled_annotations[f"{taxon}:{prot_id}"] = gt_terms
                if prot_id in prop_preds:
                    pooled_predictions[f"{taxon}:{prot_id}"] = prop_preds[prot_id]

    res = summarize_eval_lists(
        tc_all50, cc_all50, tc_gt, cc_gt, [], [], [], [],
        go, sub, pooled_annotations, pooled_predictions
    )
    res["gaef_all50_tau03"] = {
        "taxon_consistency": round(float(np.mean(tc_all50_03)), 4) if tc_all50_03 else None,
        "complex_coherence": round(float(np.mean(cc_all50_03)) / 100.0, 4) if cc_all50_03 else None,
    }
    res["gaef_gt_subset_tau03"] = {
        "taxon_consistency": round(float(np.mean(tc_gt_03)), 4) if tc_gt_03 else None,
        "complex_coherence": round(float(np.mean(cc_gt_03)) / 100.0, 4) if cc_gt_03 else None,
    }
    return res


def summarize_eval_lists(
    tc_all50, cc_all50, tc_gt, cc_gt, macro_fmax, macro_smin, macro_auc, macro_aupr,
    go, sub, pooled_annotations, pooled_predictions
):
    pooled_cafa = {}
    if pooled_annotations:
        pm = evaluate_prediction_metrics(go, sub, pooled_annotations, pooled_predictions)
        if pm is not None:
            pooled_cafa = {
                "fmax": round(float(pm["fmax"]), 4),
                "smin": round(float(pm["smin"]), 2),
                "avg_auc": round(float(pm["avg_auc"]), 4),
                "aupr": round(float(pm["aupr"]), 4),
                "wfmax": round(float(pm["wfmax"]), 4),
                "tmax": round(float(pm["tmax"]), 2),
            }
    return {
        "n_total_organisms": len(tc_all50),
        "n_gt_organisms": len(tc_gt),
        "gaef_all50": {
            "taxon_consistency": round(float(np.mean(tc_all50)), 4) if tc_all50 else None,
            "complex_coherence": round(float(np.mean(cc_all50)) / 100.0, 4) if cc_all50 else None,
        },
        "gaef_gt_subset": {
            "taxon_consistency": round(float(np.mean(tc_gt)), 4) if tc_gt else None,
            "complex_coherence": round(float(np.mean(cc_gt)) / 100.0, 4) if cc_gt else None,
        },
        "pooled_cafa": pooled_cafa,
        "macro_cafa": {
            "fmax": round(float(np.mean(macro_fmax)), 4) if macro_fmax else None,
            "smin": round(float(np.mean(macro_smin)), 2) if macro_smin else None,
            "avg_auc": round(float(np.mean(macro_auc)), 4) if macro_auc else None,
            "aupr": round(float(np.mean(macro_aupr)), 4) if macro_aupr else None,
        },
    }


def _eval_soft_single_organism_worker(task):
    taxon, fold, hpath, rpath, sub, tau, lam = task
    ctx = _WORKER_CTX
    go = ctx["go"]
    taxa_constraints_file = ctx["taxa_constraints_file"]
    homodimer_terms_gaef = ctx["homodimer_terms_gaef"]
    complex_anc_map = ctx["complex_anc_map"] if sub == "cc" else None
    all_het_complexes = ctx["all_het_complexes"] if sub == "cc" else set()
    ann_dir = ctx["ann_dir"]

    ann_file = os.path.join(ann_dir, f"annots_taxon_{taxon}.tsv")
    org_annots = load_annotations(ann_file) if os.path.exists(ann_file) else None
    gt_prot_ids = set(org_annots.keys()) if org_annots else None
    has_comp = bool(complex_anc_map)

    forb_set = ctx["naive_forb_prop"].get(taxon, set())
    tau_cut_str = "0.299999" if abs(tau - 0.3) < 1e-6 else ("0.001" if abs(tau - 0.001) < 1e-6 else None)

    union_s1 = set()
    comp_counts_s1 = defaultdict(int)
    gt_preds_s1 = {}
    above_tau_prots = defaultdict(list) if sub == "cc" else None
    best_below_tau = {} if sub == "cc" else None
    prot_complex_terms = {} if sub == "cc" else None

    with open(rpath, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            prot = parts[0]
            keep_gt = gt_prot_ids is not None and prot in gt_prot_ids
            gt_d = {} if keep_gt else None
            prot_c = set() if has_comp else None
            prot_ct = set() if has_comp else None

            for tok in parts[1:]:
                if "|" not in tok:
                    continue
                go_id, sc_str = tok.split("|", 1)
                ca = complex_anc_map.get(go_id) if has_comp else None
                is_het = (sub == "cc" and go_id in all_het_complexes)
                if not keep_gt and not is_het:
                    if tau_cut_str is not None and len(sc_str) >= 3 and sc_str[:6] <= tau_cut_str:
                        continue
                    if ca is None and (go_id in union_s1):
                        continue
                sc = float(sc_str)
                if go_id not in forb_set:
                    new_sc = sc
                else:
                    new_sc = (sc - lam) if (sc - tau) >= lam else 0.0
                if new_sc > tau:
                    union_s1.add(go_id)
                    if ca:
                        prot_c.update(ca)
                        prot_ct.add(go_id)
                if keep_gt and new_sc > 0.0:
                    gt_d[go_id] = new_sc
                if is_het and new_sc > 0.0:
                    if new_sc > tau:
                        above_tau_prots[go_id].append((prot, new_sc))
                    else:
                        prev = best_below_tau.get(go_id)
                        if prev is None or new_sc > prev[1]:
                            best_below_tau[go_id] = (prot, new_sc)

            if has_comp and prot_c:
                for c in prot_c:
                    comp_counts_s1[c] += 1
                prot_complex_terms[prot] = prot_ct
            if keep_gt:
                gt_preds_s1[prot] = gt_d

    is_sat1, c_coh1 = fast_gaef_eval_from_unions(
        union_s1, comp_counts_s1, sub, go, taxa_constraints_file, homodimer_terms_gaef
    )

    has_gt1 = False
    prop_s1 = {}
    if org_annots:
        go_set = set(go.get_namespace_terms(NAMESPACES[sub]))
        go_set.discard({"cc": "GO:0005575", "mf": "GO:0003674", "bp": "GO:0008150"}[sub])
        if any((terms & go_set) for terms in org_annots.values()):
            has_gt1 = True
            prop_s1 = propagate_scores_dict(gt_preds_s1, go, sub, gt_prot_ids)

    is_sat2, c_coh2, has_gt2, prop_s2 = None, None, False, {}
    if sub == "cc":
        union_s2 = set(union_s1)
        gt_preds_s2 = {p: dict(gd) for p, gd in gt_preds_s1.items()}
        for c_term, plist in above_tau_prots.items():
            if len(plist) == 1:
                p1, s1_val = plist[0]
                m1 = s1_val - tau
                cand = best_below_tau.get(c_term)
                if cand is not None:
                    p2, s2_val = cand
                    d2 = tau - s2_val
                else:
                    p2, d2 = None, float("inf")
                if min(m1, d2) < lam:
                    if m1 <= d2:
                        if p1 in prot_complex_terms:
                            prot_complex_terms[p1].discard(c_term)
                        if p1 in gt_preds_s2 and c_term in gt_preds_s2[p1]:
                            gt_preds_s2[p1][c_term] = min(gt_preds_s2[p1][c_term], tau)
                    elif p2 is not None:
                        union_s2.add(c_term)
                        prot_complex_terms.setdefault(p2, set()).add(c_term)
                        if p2 in gt_preds_s2:
                            gt_preds_s2[p2][c_term] = tau + 0.001

        comp_counts_s2 = defaultdict(int)
        for p, ct_set in prot_complex_terms.items():
            pc = set()
            for g in ct_set:
                ca = complex_anc_map.get(g)
                if ca:
                    pc.update(ca)
            for c in pc:
                comp_counts_s2[c] += 1

        is_sat2, c_coh2 = fast_gaef_eval_from_unions(
            union_s2, comp_counts_s2, sub, go, taxa_constraints_file, homodimer_terms_gaef
        )
        if has_gt1:
            has_gt2 = True
            prop_s2 = propagate_scores_dict(gt_preds_s2, go, sub, gt_prot_ids)

    return (
        taxon,
        is_sat1, c_coh1, has_gt1, org_annots if has_gt1 else None, prop_s1 if has_gt1 else None,
        is_sat2, c_coh2, has_gt2, prop_s2 if has_gt2 else None,
    )


def evaluate_soft_regularizer_on_timeset(
    raw_dir: str,
    hard_stage1_dir: str,
    sub: str,
    tau: float,
    lam: float,
    timeset_taxons: set,
    ann_dir: str,
    go: Ontology,
    go_is_a: Ontology,
    taxa_constraints_file: str,
    complex_child_terms: set,
    homodimer_terms_gaef: set,
    homodimer_terms_s2: set = None,
    go_hierarchy_cc: dict = None,
    pool=None,
):
    """
    Evaluate Soft Regularizer (SBR) at penalty lambda for both Stage 1 (soft_s1)
    and (when sub == 'cc') Stage 2 (soft_s2) in parallel across organisms.
    """
    tc_s1_all, cc_s1_all, tc_s1_gt, cc_s1_gt = [], [], [], []
    pool_ann_s1, pool_pred_s1 = {}, {}

    tc_s2_all, cc_s2_all, tc_s2_gt, cc_s2_gt = [], [], [], []
    pool_ann_s2, pool_pred_s2 = {}, {}

    tasks = []
    for rpath in sorted(glob.glob(os.path.join(raw_dir, "predictions_fold_*_taxon_*.tsv"))):
        m = FN_PRED.search(os.path.basename(rpath))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        if taxon not in timeset_taxons:
            continue
        tasks.append((taxon, fold, None, rpath, sub, tau, lam))

    results = list(pool.map(_eval_soft_single_organism_worker, tasks)) if pool is not None else [
        _eval_soft_single_organism_worker(t) for t in tasks
    ]

    for (
        taxon,
        is_sat1, c_coh1, has_gt1, org_annots, prop_s1,
        is_sat2, c_coh2, has_gt2, prop_s2,
    ) in results:
        tc_s1_all.append(is_sat1)
        if c_coh1 is not None:
            cc_s1_all.append(c_coh1)
        if has_gt1:
            tc_s1_gt.append(is_sat1)
            if c_coh1 is not None:
                cc_s1_gt.append(c_coh1)
            for prot_id, gt_terms in org_annots.items():
                pool_ann_s1[f"{taxon}:{prot_id}"] = gt_terms
                if prot_id in prop_s1:
                    pool_pred_s1[f"{taxon}:{prot_id}"] = prop_s1[prot_id]

        if sub == "cc":
            tc_s2_all.append(is_sat2)
            if c_coh2 is not None:
                cc_s2_all.append(c_coh2)
            if has_gt2:
                tc_s2_gt.append(is_sat2)
                if c_coh2 is not None:
                    cc_s2_gt.append(c_coh2)
                for prot_id, gt_terms in org_annots.items():
                    pool_ann_s2[f"{taxon}:{prot_id}"] = gt_terms
                    if prot_id in prop_s2:
                        pool_pred_s2[f"{taxon}:{prot_id}"] = prop_s2[prot_id]

    res_s1 = summarize_eval_lists(
        tc_s1_all, cc_s1_all, tc_s1_gt, cc_s1_gt, [], [], [], [],
        go, sub, pool_ann_s1, pool_pred_s1
    )
    res_s2 = (
        summarize_eval_lists(
            tc_s2_all, cc_s2_all, tc_s2_gt, cc_s2_gt, [], [], [], [],
            go, sub, pool_ann_s2, pool_pred_s2
        )
        if sub == "cc"
        else None
    )
    return res_s1, res_s2


def run_stage2_on_dir(
    stage1_dir: str,
    stage2_out_dir: str,
    complexes_file: str,
    go_hierarchy_cc_file: str,
    tau: float,
    demotion_only: bool,
    target_taxons: set = None,
):
    os.makedirs(stage2_out_dir, exist_ok=True)
    total_flips = 0
    total_annot = 0
    n_active_orgs = 0
    n_total_orgs = 0

    for fpath in sorted(glob.glob(os.path.join(stage1_dir, "optimized_fold_*_taxon_*.tsv"))):
        m = FN_PRED.search(os.path.basename(fpath))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        if target_taxons is not None and taxon not in target_taxons:
            continue
        n_total_orgs += 1
        out_path = os.path.join(stage2_out_dir, os.path.basename(fpath))
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            continue
        flips, n_ann = complex_coherence_adjust(
            predictions_file=fpath,
            complexes_file=complexes_file,
            go_hierarchy_file=go_hierarchy_cc_file,
            output_file=out_path,
            threshold=tau,
            optimized=True,
            top_k=50,
            demotion_only=demotion_only,
        )
        if n_ann > 0:
            n_active_orgs += 1
        total_flips += flips
        total_annot += n_ann

    return {
        "n_total_orgs": n_total_orgs,
        "n_active_incoherent_orgs": n_active_orgs,
        "total_flips": total_flips,
        "total_annotated_predictions_in_active": total_annot,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_dir", required=True)
    ap.add_argument("--code_dir", required=True)
    ap.add_argument("--gaef_dir", required=True)
    ap.add_argument("--timeset_taxons_file", required=True)
    ap.add_argument("--stage1_root", required=True)
    ap.add_argument("--funtaxis_root", required=True)
    ap.add_argument("--naive_root", required=True)
    ap.add_argument("--predictor", required=True, choices=["seq_sim_timeset", "deepgo-se", "sprof-go"])
    ap.add_argument("--sub", default="all", choices=["all", "cc", "mf", "bp"])
    ap.add_argument("--output_json", required=True)
    args = ap.parse_args()

    with open(args.timeset_taxons_file) as fh:
        timeset_taxons = {line.strip() for line in fh if line.strip()}

    go_obo_file = os.path.join(args.code_dir, "data", "go-basic.obo")
    go = Ontology(go_obo_file)
    go_is_a = Ontology(go_obo_file, with_rels=False)
    term_to_children, _, _ = coherence.parse_go_ontology(go_obo_file, include_part_of=False)
    complex_child_terms = coherence.get_all_child_terms("GO:0032991", term_to_children)
    complex_child_terms.add("GO:0032991")

    cdir = os.path.join(args.code_dir, "data", "constraints")
    taxa_constraints_file = os.path.join(cdir, "taxon_constraints.tsv")
    homodimer_terms_gaef = parse_homodimer_terms_file(os.path.join(cdir, "protein_complexes.tsv"))

    ann_dir = os.path.join(args.data_dir, "swissprot_proteomes_folds", "annotations-go-basic")
    complexes_file = os.path.join(args.code_dir, "complex_coherence", "protein_complexes.tsv")
    go_hierarchy_cc_file = os.path.join(args.code_dir, "data", "go_hierarchy_cc.tsv")
    homodimer_terms_s2 = load_homodimer_terms(complexes_file)
    go_hierarchy_cc = load_cc_go_hierarchy(go_hierarchy_cc_file)

    # Precompute naive_noprop and naive_prop forbidden sets per timeset taxon (matching naive_filter.py)
    in_taxon, never_in = load_constraints(os.path.join(args.code_dir, "data", "go_taxon_constraints_extracted_obo.tsv"))
    go_parents = load_full_go_hierarchy(os.path.join(args.code_dir, "data", "go_hierarchy.tsv"))
    go_children = defaultdict(set)
    for c, ps in go_parents.items():
        for p in ps:
            go_children[p].add(c)
    tax_isa, _, _ = load_taxon_hierarchy(
        os.path.join(args.code_dir, "data", "taxon_hierarchy.tsv"),
        os.path.join(args.code_dir, "data", "ncbitaxon_hierarchy.tsv"),
    )
    naive_forb_noprop = {}
    naive_forb_prop = {}
    for tnum in timeset_taxons:
        lineage = transitive({f"NCBITaxon_{tnum}"}, tax_isa)
        violated = set()
        for g, taxa in never_in.items():
            if any(t in lineage for t in taxa):
                violated.add(g)
        for g, taxa in in_taxon.items():
            if taxa and not any(t in lineage for t in taxa):
                violated.add(g)
        naive_forb_noprop[tnum] = set(violated)
        naive_forb_prop[tnum] = transitive(violated, go_children)

    default_taus = {
        "seq_sim_timeset": 0.001,
        "deepgo-se": 0.3,
        "sprof-go": 0.001,
    }
    tau = default_taus[args.predictor]

    complex_anc_map = {}
    for g in go.get_namespace_terms(NAMESPACES["cc"]):
        ca = _get_complex_ancestors(g, go_is_a, complex_child_terms)
        if ca:
            complex_anc_map[g] = ca
    all_het_complexes = get_heteromeric_complexes(homodimer_terms_s2, set(go_hierarchy_cc.keys()), go_hierarchy_cc)

    import multiprocessing as mp
    n_workers = min(16, os.cpu_count() or 16)
    worker_ctx = {
        "go": go,
        "taxa_constraints_file": taxa_constraints_file,
        "homodimer_terms_gaef": homodimer_terms_gaef,
        "complex_anc_map": complex_anc_map,
        "all_het_complexes": all_het_complexes,
        "ann_dir": ann_dir,
        "naive_forb_noprop": naive_forb_noprop,
        "naive_forb_prop": naive_forb_prop,
    }

    if os.path.exists(args.output_json):
        with open(args.output_json) as fh:
            out_summary = json.load(fh)
        out_summary["predictor"] = args.predictor
        out_summary["solver_tau"] = tau
        out_summary.setdefault("subontologies", {})
    else:
        out_summary = {"predictor": args.predictor, "solver_tau": tau, "subontologies": {}}
    os.makedirs(os.path.dirname(args.output_json), exist_ok=True)

    target_subs = ["cc", "mf", "bp"] if args.sub == "all" else [args.sub]
    with mp.get_context("fork").Pool(processes=n_workers, initializer=_init_worker, initargs=(worker_ctx,)) as pool:
        for sub in target_subs:
            out_summary["subontologies"].setdefault(sub, {})
            raw_dir = os.path.join(
                args.data_dir, "swissprot_proteomes_folds", f"{args.predictor}_results", sub, "predictions"
            )
            s1_with_dir = os.path.join(args.stage1_root, args.predictor, sub, "optimized")
            if not os.path.exists(s1_with_dir):
                s1_with_dir = os.path.join(
                    args.data_dir, "swissprot_proteomes_folds", f"{args.predictor}_results", sub, "optimized"
                )
            s1_notax_dir = os.path.join(args.stage1_root, args.predictor, sub, "optimized_without_taxons")
            if not os.path.exists(s1_notax_dir):
                s1_notax_dir = os.path.join(
                    args.data_dir, "swissprot_proteomes_folds", f"{args.predictor}_results", sub, "optimized_without_taxons"
                )
            ft_dir = os.path.join(args.funtaxis_root, args.predictor, sub, "funtaxis_lite")

            eval_specs = [
                ("funtaxis_lite_prop", ft_dir, None),
                ("stage1_with_taxons", s1_with_dir, None),
                ("stage1_no_taxons", s1_notax_dir, None),
                ("before", raw_dir, None),
                ("naive_noprop", raw_dir, "noprop"),
                ("naive_prop", raw_dir, "prop"),
            ]

            if sub == "cc":
                orig_s2_with = os.path.join(
                    args.data_dir, "swissprot_proteomes_folds", f"{args.predictor}_results", "cc", "optimized_complex_coherence"
                )
                orig_s2_notax = os.path.join(
                    args.data_dir, "swissprot_proteomes_folds", f"{args.predictor}_results", "cc", "optimized_without_taxons_complex_coherence"
                )
                if os.path.exists(orig_s2_with):
                    eval_specs.insert(0, ("existing_stage2_standard_withtax", orig_s2_with, None))
                if os.path.exists(orig_s2_notax):
                    eval_specs.insert(1, ("existing_stage2_standard_notax", orig_s2_notax, None))

                for s1_label, s1_d in [("withtax", s1_with_dir), ("notax", s1_notax_dir)]:
                    if not os.path.exists(s1_d):
                        continue
                    s2_dem_dir = os.path.join(args.stage1_root, args.predictor, "cc", f"stage2_demotion_only_{s1_label}")
                    if f"stage2_stats_demotion_only_{s1_label}" not in out_summary["subontologies"]["cc"]:
                        stats_dem = run_stage2_on_dir(
                            s1_d, s2_dem_dir, complexes_file, go_hierarchy_cc_file, tau, True, timeset_taxons
                        )
                        out_summary["subontologies"]["cc"][f"stage2_stats_demotion_only_{s1_label}"] = stats_dem
                    eval_specs.insert(0, (f"stage2_demotion_only_{s1_label}", s2_dem_dir, None))

                    if abs(tau - 0.3) > 1e-6:
                        s2_dem03_dir = os.path.join(args.stage1_root, args.predictor, "cc", f"stage2_demotion_only_tau03_{s1_label}")
                        if f"stage2_stats_demotion_only_tau03_{s1_label}" not in out_summary["subontologies"]["cc"]:
                            stats_dem03 = run_stage2_on_dir(
                                s1_d, s2_dem03_dir, complexes_file, go_hierarchy_cc_file, 0.3, True, timeset_taxons
                            )
                            out_summary["subontologies"]["cc"][f"stage2_stats_demotion_only_tau03_{s1_label}"] = stats_dem03
                        eval_specs.insert(0, (f"stage2_demotion_only_tau03_{s1_label}", s2_dem03_dir, None))

                    if s1_label == "withtax":
                        s2_std_dir = os.path.join(args.stage1_root, args.predictor, "cc", f"stage2_standard_{s1_label}")
                        if f"stage2_stats_standard_{s1_label}" not in out_summary["subontologies"]["cc"]:
                            stats_std = run_stage2_on_dir(
                                s1_d, s2_std_dir, complexes_file, go_hierarchy_cc_file, tau, False, timeset_taxons
                            )
                            out_summary["subontologies"]["cc"][f"stage2_stats_standard_{s1_label}"] = stats_std
                        eval_specs.insert(0, (f"stage2_standard_{s1_label}", s2_std_dir, None))

            for label, dpath, forb_map in eval_specs:
                if label in out_summary["subontologies"][sub]:
                    continue
                if not os.path.exists(dpath):
                    continue
                res_tau = evaluate_dir_on_timeset(
                    dpath,
                    sub,
                    tau,
                    timeset_taxons,
                    ann_dir,
                    go,
                    go_is_a,
                    taxa_constraints_file,
                    complex_child_terms,
                    homodimer_terms_gaef,
                    forbidden_by_taxon=forb_map,
                    pool=pool,
                )
                out_summary["subontologies"][sub][label] = res_tau
                print(f"[{args.predictor} {sub} {label}] {json.dumps(res_tau)}", flush=True)
                with open(args.output_json, "w") as fh:
                    json.dump(out_summary, fh, indent=2)

            # Evaluate Soft Regularizer (SBR) at lambda=0.2 and lambda=0.5 in-memory
            for lam in [0.2, 0.5]:
                if f"soft_s1_lam{lam}" in out_summary["subontologies"][sub] and (
                    sub != "cc" or f"soft_s2_lam{lam}" in out_summary["subontologies"][sub]
                ):
                    continue
                res_s1, res_s2 = evaluate_soft_regularizer_on_timeset(
                    raw_dir,
                    s1_with_dir,
                    sub,
                    tau,
                    lam,
                    timeset_taxons,
                    ann_dir,
                    go,
                    go_is_a,
                    taxa_constraints_file,
                    complex_child_terms,
                    homodimer_terms_gaef,
                    homodimer_terms_s2=homodimer_terms_s2 if sub == "cc" else None,
                    go_hierarchy_cc=go_hierarchy_cc if sub == "cc" else None,
                    pool=pool,
                )
                out_summary["subontologies"][sub][f"soft_s1_lam{lam}"] = res_s1
                print(f"[{args.predictor} {sub} soft_s1_lam{lam}] {json.dumps(res_s1)}", flush=True)
                if res_s2 is not None:
                    out_summary["subontologies"][sub][f"soft_s2_lam{lam}"] = res_s2
                    print(f"[{args.predictor} {sub} soft_s2_lam{lam}] {json.dumps(res_s2)}", flush=True)
                with open(args.output_json, "w") as fh:
                    json.dump(out_summary, fh, indent=2)


if __name__ == "__main__":
    main()
