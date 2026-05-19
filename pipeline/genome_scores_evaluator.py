"""
Given a taxon scores, evaluate metrics. Based on GAEF framework.

Metrics:
1. Completeness
2. Coherence
3. Consistency

4. F-max
5. Precision
6. Recall
7. S-min
8. AUC
9. AUPR
10. AVGIC

"""

import sys
import os

# Make the repo root importable so `gaef_patches` and `pipeline` resolve, and
# the upstream GAEF clone importable for `GAEF.completeness` / `GAEF.utils`.
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
# Upstream GAEF clone is expected as a sibling of the repo (or at ${GAEF_DIR}).
_GAEF_PARENT = os.environ.get('GAEF_DIR') and os.path.dirname(os.path.abspath(os.environ['GAEF_DIR']))
for _p in (
    os.path.abspath(os.path.join(_REPO_ROOT, '..')),
    _GAEF_PARENT,
):
    if _p and _p not in sys.path:
        sys.path.append(_p)

import csv
import json
import argparse
from collections import defaultdict
from typing import Set, DefaultDict, Optional
from pathlib import Path
import numpy as np
import pandas as pd

# Upstream-clean GAEF modules.
from GAEF.completeness import read_essential_terms, count_essential_terms
from GAEF.utils import Ontology, FUNC_DICT, NAMESPACES

from gaef_patches.taxon_consistency import TaxonHierarchy, parse_owl_files
from gaef_patches.complex_classifier import parse_homodimer_terms_file
from gaef_patches import coherence

# Vendored constraint files live under <repo>/data/constraints/.
_DEFAULT_CONSTRAINTS_DIR = os.path.join(_REPO_ROOT, 'data', 'constraints')

from sklearn.metrics import roc_curve, auc
import math


def evaluate_annotations_optimized(go, real_annots, pred_annots, ic_cache=None):
    """
    Optimized version of evaluate_annotations with IC value caching.

    Args:
       go (utils.Ontology): Ontology class instance with go.obo
       real_annots (list): List of sets of real GO classes
       pred_annots (list): List of sets of predicted GO classes
       ic_cache (dict): Optional pre-computed cache of {term: (ic, norm_ic)}

    Returns:
       Tuple of (f, p, r, s, ru, mi, fps, fns, avg_ic, wf)
    """
    if ic_cache is None:
        ic_cache = {}

    total = 0
    p = 0.0
    r = 0.0
    wp = 0.0
    wr = 0.0
    p_total = 0
    ru = 0.0
    mi = 0.0
    avg_ic = 0.0
    fps = []
    fns = []

    for i in range(len(real_annots)):
        if len(real_annots[i]) == 0:
            continue

        tp = real_annots[i] & pred_annots[i]  # Set intersection (faster than set().intersection())
        fp = pred_annots[i] - tp
        fn = real_annots[i] - tp

        tpic = 0.0
        for go_id in tp:
            if go_id not in ic_cache:
                ic_cache[go_id] = (go.get_ic(go_id), go.get_norm_ic(go_id))
            ic_val, norm_ic_val = ic_cache[go_id]
            tpic += norm_ic_val
            avg_ic += ic_val

        fpic = 0.0
        for go_id in fp:
            if go_id not in ic_cache:
                ic_cache[go_id] = (go.get_ic(go_id), go.get_norm_ic(go_id))
            ic_val, norm_ic_val = ic_cache[go_id]
            fpic += norm_ic_val
            mi += ic_val

        fnic = 0.0
        for go_id in fn:
            if go_id not in ic_cache:
                ic_cache[go_id] = (go.get_ic(go_id), go.get_norm_ic(go_id))
            ic_val, norm_ic_val = ic_cache[go_id]
            fnic += norm_ic_val
            ru += ic_val

        fps.append(fp)
        fns.append(fn)
        tpn = len(tp)
        fpn = len(fp)
        fnn = len(fn)
        total += 1
        recall = tpn / (1.0 * (tpn + fnn))
        r += recall
        wrecall = tpic / (tpic + fnic) if (tpic + fnic) > 0 else 0.0
        wr += wrecall
        if len(pred_annots[i]) > 0:
            p_total += 1
            precision = tpn / (1.0 * (tpn + fpn))
            p += precision
            if tpic + fpic > 0:
                wp += tpic / (tpic + fpic)

    if total == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, fps, fns, 0.0, 0.0

    avg_ic = (avg_ic + mi) / total
    ru /= total
    mi /= total
    r /= total
    wr /= total
    if p_total > 0:
        p /= p_total
        wp /= p_total
    f = 0.0
    wf = 0.0
    if p + r > 0:
        f = 2 * p * r / (p + r)
    if wp + wr > 0:
        wf = 2 * wp * wr / (wp + wr)
    s = math.sqrt(ru * ru + mi * mi)
    return f, p, r, s, ru, mi, fps, fns, avg_ic, wf


def compute_roc_auc(labels, preds):
    """Compute ROC AUC for a single term."""
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    return auc(fpr, tpr)


def compute_metrics_optimized(test_df, go, terms_dict, terms, ont, eval_preds, labels_matrix=None):
    """
    Optimized version of compute_metrics with:
    - Pre-computed labels_matrix support
    - Set intersection instead of map/lambda/filter
    - IC value caching

    Args:
        test_df: DataFrame with proteins, preds, prop_annotations, exp_annotations
        go: Ontology object
        terms_dict: Dict mapping term -> index
        terms: List of terms
        ont: Ontology namespace ('cc', 'bp', 'mf')
        eval_preds: numpy array of prediction scores (num_proteins x num_terms)
        labels_matrix: Optional pre-computed binary labels matrix (num_proteins x num_terms)

    Returns:
        Tuple of (fmax, smin, tmax, wfmax, wtmax, avg_auc, aupr, avgic, fmax_spec_match, precision_at_tmax, recall_at_tmax, tp_at_tmax, fp_at_tmax, fn_at_tmax, sensitivity_at_tmax, specificity_at_tmax)
    """
    num_proteins = len(test_df)
    num_terms = len(terms_dict)

    # Build labels matrix if not provided
    if labels_matrix is None:
        labels_matrix = np.zeros((num_proteins, num_terms), dtype=np.float32)
        for i, row in enumerate(test_df.itertuples()):
            for go_id in row.prop_annotations:
                if go_id in terms_dict:
                    labels_matrix[i, terms_dict[go_id]] = 1

    # Compute per-term AUC
    total_n = 0
    total_sum = 0.0
    for i in range(num_terms):
        pos_n = np.sum(labels_matrix[:, i])
        if pos_n > 0 and pos_n < num_proteins:
            total_n += 1
            roc_auc = compute_roc_auc(labels_matrix[:, i], eval_preds[:, i])
            total_sum += roc_auc

    avg_auc = total_sum / total_n if total_n > 0 else 0.0

    # Prepare labels using set intersection (faster than map/lambda/filter)
    go_set = go.get_namespace_terms(NAMESPACES[ont])
    root_term = FUNC_DICT[ont]
    if root_term in go_set:
        go_set.remove(root_term)

    # Use set intersection instead of map/lambda/filter
    labels = [annots & go_set for annots in test_df['prop_annotations'].values]
    spec_labels = [annots & go_set for annots in test_df['exp_annotations'].values]

    # Initialize tracking variables
    fmax = 0.0
    tmax = 0.0
    wfmax = 0.0
    wtmax = 0.0
    avgic = 0.0
    smin = 1000000.0
    fmax_spec_match = 0
    precisions = []
    recalls = []
    precision_at_tmax = 0.0
    recall_at_tmax = 0.0
    tp_at_tmax = 0
    fp_at_tmax = 0
    fn_at_tmax = 0

    # IC cache for evaluate_annotations_optimized
    ic_cache = {}

    # Threshold sweep
    for t in range(0, 101):
        threshold = t / 100.0

        # Vectorized threshold application
        above_threshold = eval_preds >= threshold

        # Build prediction sets using set intersection with go_set
        preds = []
        for i in range(num_proteins):
            above_indices = np.where(above_threshold[i])[0]
            annots = {terms[j] for j in above_indices}
            preds.append(annots & go_set)  # Set intersection instead of filter

        # Evaluate at this threshold
        fscore, prec, rec, s, ru, mi, fps, fns, avg_ic, wf = evaluate_annotations_optimized(
            go, labels, preds, ic_cache)

        # Count spec matches
        spec_match = sum(len(spec_labels[i] & preds[i]) for i in range(num_proteins))

        # Aggregate TP, FP, FN (fps/fns align with labels indices for non-empty real_annots)
        tp_sum = sum(len(labels[i] & preds[i]) for i in range(num_proteins) if len(labels[i]) > 0)
        fp_sum = sum(len(fp) for fp in fps)
        fn_sum = sum(len(fn) for fn in fns)

        precisions.append(prec)
        recalls.append(rec)

        if fmax < fscore:
            fmax = fscore
            tmax = threshold
            avgic = avg_ic
            fmax_spec_match = spec_match
            precision_at_tmax = prec
            recall_at_tmax = rec
            tp_at_tmax = tp_sum
            fp_at_tmax = fp_sum
            fn_at_tmax = fn_sum
        if wfmax < wf:
            wfmax = wf
            wtmax = threshold
        if smin > s:
            smin = s

    # Compute AUPR
    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    aupr = np.trapz(precisions, recalls)

    # Term-centric sensitivity and specificity at F-max threshold
    above_at_tmax = eval_preds >= tmax
    sens_list, spec_list = [], []
    for i in range(num_terms):
        labels_i = labels_matrix[:, i]
        preds_i = above_at_tmax[:, i]
        tp = np.sum((labels_i == 1) & preds_i)
        fn = np.sum((labels_i == 1) & ~preds_i)
        fp = np.sum((labels_i == 0) & preds_i)
        tn = np.sum((labels_i == 0) & ~preds_i)
        pos_n = tp + fn
        neg_n = tn + fp
        if pos_n > 0:
            sens_list.append(tp / pos_n)
        if neg_n > 0:
            spec_list.append(tn / neg_n)
    sensitivity_at_tmax = float(np.mean(sens_list)) if sens_list else 0.0
    specificity_at_tmax = float(np.mean(spec_list)) if spec_list else 0.0

    return fmax, smin, tmax, wfmax, wtmax, avg_auc, aupr, avgic, fmax_spec_match, precision_at_tmax, recall_at_tmax, tp_at_tmax, fp_at_tmax, fn_at_tmax, sensitivity_at_tmax, specificity_at_tmax


def get_ancestors_with_ontology(protein_go_terms, go):
    """
    Expands GO terms to include their ancestors for each protein.
    Optimized version with ancestor caching for large datasets.
    
    Parameters:
    - protein_go_terms (dict): {protein_id: set(GO_terms)}
    - go: Ontology object (reused)
    
    Returns:
    - dict: {protein_id: set(GO_terms + ancestors)}
    """
    result = {}
    # Cache ancestors to avoid recomputation across all proteins
    ancestor_cache = {}
    
    for protein_id, terms in protein_go_terms.items():
        expanded = set()
        for term in terms:
            # Use cache to avoid redundant get_ancestors calls
            if term not in ancestor_cache:
                ancestor_cache[term] = go.get_ancestors(term)
            expanded.update(ancestor_cache[term])
        result[protein_id] = expanded
    return result, ancestor_cache

def get_specific_with_ontology(protein_go_terms, ontology, ancestor_cache=None):
    """
    Get specific (non-redundant) GO terms for each protein.
    Optimized version with ancestor caching for large datasets.

    Parameters:
    - protein_go_terms (dict): {protein_id: set(GO_terms)}
    - ontology: Ontology object
    - ancestor_cache: Optional pre-computed cache of {term: ancestors} from get_ancestors_with_ontology

    Returns:
    - dict: {protein_id: set(specific GO_terms)}
    """
    if ancestor_cache is None:
        ancestor_cache = {}

    result = {}
    for protein_id, terms in protein_go_terms.items():
        # Collect all ancestors for this protein's terms using cache
        all_ancestors = set()
        for term in terms:
            # Use cache to avoid redundant get_ancestors calls
            if term not in ancestor_cache:
                ancestor_cache[term] = ontology.get_ancestors(term)
            ancestors = ancestor_cache[term]
            # Add ancestors excluding the term itself
            all_ancestors.update(ancestors)
            all_ancestors.discard(term)
        # Single set subtraction - keep only terms that aren't ancestors of other terms
        result[protein_id] = terms - all_ancestors
    return result

def load_predictions(predictions_file, threshold=0.5):
    """
    Load predictions from TSV file, optimized to parse strings only once.
    
    Args:
        predictions_file: Path to TSV file
        threshold: Score threshold for filtering
    
    Returns:
        Dictionary mapping protein_id -> set(go_terms)
    """
    protein_go_terms = {}
    with open(predictions_file, 'r') as f:
        tsv_reader = csv.reader(f, delimiter='\t')
        for row in tsv_reader:
            if row:
                protein_id = row[0]
                protein_go_terms[protein_id] = {go_term_score.split('|')[0] for go_term_score in row[1:] if float(go_term_score.split('|')[1]) > threshold}
    return protein_go_terms


def check_consistency_fast(protein_go_terms, taxa_constraints_file, owl_dir=None):
    """
    Check taxonomic consistency of protein-GO term annotations using a fast implementation.

    Args:
        protein_go_terms: Dictionary mapping protein_id -> set of GO terms
        taxa_constraints_file: Path to TSV file with taxon constraints
        owl_dir: Directory containing `ncbitaxon_with_disjointness.owl` and
                 `go-taxon-groupings.owl`. Defaults to the directory containing
                 `taxa_constraints_file`.

    Returns:
        True if consistent (satisfiable), False otherwise
    """
    # === Process constraints ===
    constraints = pd.read_csv(taxa_constraints_file, sep="\t", dtype=str)
    only_map = defaultdict(set)
    never_map = defaultdict(set)

    for _, row in constraints.iterrows():
        go_id = row["GO_ID"]
        lineage = row["Taxon_ID"]
        if row["Constraint_Type"] == "only_in_taxon":
            only_map[go_id].add(lineage)  # expected to be normalized: NCBITaxon_<taxon_id>
        elif row["Constraint_Type"] == "never_in_taxon":
            never_map[go_id].add(lineage)  # expected to be normalized: NCBITaxon_<taxon_id>

    # Build GenomeAnnotation objects in memory
    annotations = defaultdict(dict)
    for protein_id, go_terms in protein_go_terms.items():
        for go_id in go_terms:
            if go_id in annotations:
                continue
            only_taxa = only_map.get(go_id, set())
            never_taxa = never_map.get(go_id, set())
            if only_taxa or never_taxa:
                annotations[go_id] = {
                    "never_in_taxon": never_taxa,
                    "only_in_taxon": only_taxa
                }
    
    # Early exit if no constraints
    if not annotations:
        return True
    
    # === Load taxon hierarchy ===
    owl_dir_path = Path(owl_dir) if owl_dir is not None else Path(taxa_constraints_file).parent
    taxon_file = owl_dir_path / "ncbitaxon_with_disjointness.owl"
    go_taxon_file = owl_dir_path / "go-taxon-groupings.owl"
    
    # Parse ontologies
    subclass_of, disjoint_with, union_members = parse_owl_files(
        str(taxon_file), str(go_taxon_file)
    )
    
    # Build hierarchy
    hierarchy = TaxonHierarchy(subclass_of, disjoint_with, union_members)
    
    # === Fast satisfiability check ===
    return _is_satisfiable_fast(annotations, hierarchy)


def _is_satisfiable_fast(annotations: DefaultDict[str, DefaultDict[str, Set[str]]], hierarchy: TaxonHierarchy) -> bool:
    """
    Fast satisfiability check without explanation generation.

    Returns:
        True if satisfiable, False otherwise
    """
    # Collect constraints - use direct set operations
    only_set: Set[str] = set()
    never_set: Set[str] = set()

    for go_id, constraints in annotations.items():
        for taxon in constraints["only_in_taxon"]:
            only_set.add(taxon.replace("NCBITaxon_", ""))
        for taxon in constraints["never_in_taxon"]:
            never_set.add(taxon.replace("NCBITaxon_", ""))

    # Early exit if no constraints
    if not only_set and not never_set:
        return True

    # Check 1: Direct conflict - same taxon in both
    if only_set & never_set:
        return False

    # Check 2: Mutual compatibility of "only_in" constraints
    only_list = list(only_set)
    for i, taxon_a in enumerate(only_list):
        full_a = f"NCBITaxon_{taxon_a}"
        for taxon_b in only_list[i+1:]:
            full_b = f"NCBITaxon_{taxon_b}"
            if hierarchy.are_disjoint(full_a, full_b)[0]:
                return False

    # Check 3: "never_in" vs "only_in" conflicts
    for only_taxon in only_set:
        full_only = f"NCBITaxon_{only_taxon}"
        only_ancestors = hierarchy.get_ancestors(full_only)

        for never_taxon in never_set:
            full_never = f"NCBITaxon_{never_taxon}"

            # If only_taxon is a subclass of never_taxon, conflict
            if full_never in only_ancestors:
                return False

    # Check 4: "only" disjoint with "never_neg"
    for only_taxon in only_set:
        full_only = f"NCBITaxon_{only_taxon}"

        for never_taxon in never_set:
            neg_class = f"NCBITaxon_{never_taxon}_neg"
            if hierarchy.are_disjoint(full_only, neg_class)[0]:
                return False

    return True

def compute_ic_depth_breadth_fast(protein_go_terms_specific: DefaultDict[str, Set[str]], ontology: Ontology):
    """
    Compute IC depth and breadth for a given set of protein-GO terms.
    
    - ic_depth: average IC per annotation
    - ic_breadth: total sum of IC values
    - normalized_ic_breadth: average IC per protein
    """

    propagated_terms_dict, _ = get_ancestors_with_ontology(protein_go_terms_specific, ontology)
    # Convert dict to list of sets for calculate_ic (expects iterable of iterables)
    propagated_terms_list = list(propagated_terms_dict.values())
    ontology.calculate_ic(propagated_terms_list)
    
    # Calculate metrics from original specific terms (not propagated)
    total_ic = 0.0
    total_annotations = 0
    num_proteins = 0
    
    for protein_id, specific_terms in protein_go_terms_specific.items():
        if specific_terms:
            for term in specific_terms:
                ic_val = ontology.get_norm_ic(term)
                total_ic += ic_val
                total_annotations += 1
            num_proteins += 1
    
    # Calculate metrics matching information_content.py formulas
    ic_depth = total_ic / total_annotations if total_annotations else 0.0
    ic_breadth = total_ic
    normalized_ic_breadth = total_ic / num_proteins if num_proteins else 0.0
    
    return ic_depth, ic_breadth, normalized_ic_breadth


def GAEF_evaluation_from_terms(assembly_name, protein_go_terms, GAEF_dir, go_file, subontology=None, output_file=None, constraints_dir=None):
    """
    Compute GAEF metrics from pre-loaded protein GO terms.
    Works with data from either load_predictions() or load_annotations().

    Args:
        assembly_name: Name of the assembly/genome
        protein_go_terms: dict[str, set[str]] mapping protein_id -> set of GO terms
        GAEF_dir: Kept for backward compatibility — no longer used to locate
                  constraint files (those are vendored under <repo>/data/constraints/).
        go_file: Path to GO OBO file
        subontology: Optional filter (cc/bp/mf, default: all)
        output_file: Optional path for JSON output (default: assembly_name + "_report.json")
        constraints_dir: Override directory holding the vendored constraint files.
                         Defaults to <repo>/data/constraints/.

    Returns:
        Dictionary with GAEF metrics (JSON-serializable)
    """
    cdir = constraints_dir if constraints_dir is not None else _DEFAULT_CONSTRAINTS_DIR
    term_file = os.path.join(cdir, "essential_terms.tsv")
    has_part_file = os.path.join(cdir, "has_part_relations.txt")
    ec2go_file = os.path.join(cdir, "ec2go_v2025-03-16")
    pathway_file = os.path.join(cdir, "metacyc_GO_v2025-03-16_with_EC.tsv")
    ontology_file = go_file
    taxa_constraints_file = os.path.join(cdir, "taxon_constraints.tsv")
    MACROMOLECULAR_COMPLEX = "GO:0032991"
    HOMODIMERIZATION = "GO:0042803"
    homodimer_terms_file = os.path.join(cdir, "protein_complexes.tsv")
    homodimer_terms = parse_homodimer_terms_file(homodimer_terms_file)

    # Load ontology once and reuse
    go = Ontology(ontology_file)

    # Expand terms to ancestors and get specific terms
    protein_go_terms_ancestors, ancestor_cache = get_ancestors_with_ontology(protein_go_terms, go)
    protein_go_terms_specific = get_specific_with_ontology(protein_go_terms, go, ancestor_cache)

    ### COMPLETENESS ###
    # Essential terms
    if subontology is None or subontology == 'bp':
        core_terms, periph_terms = read_essential_terms(term_file)
        core_ids = {e['term'] for e in core_terms}
        periph_ids = {e['term'] for e in periph_terms}

        core_count = count_essential_terms(protein_go_terms_ancestors, core_ids)
        periph_count = count_essential_terms(protein_go_terms_ancestors, periph_ids)
        essential_percentage = (sum(core_count.values())/len(core_ids)) * 100 if core_ids else 0

        # Tables grouping
        go_core   = {'Core': core_terms}
        go_periph = defaultdict(list)
        for e in periph_terms:
            go_periph[e['category']].append(e)

        # Found terms set
        found_terms = {term for term, pres in core_count.items() if pres} | {term for term, pres in periph_count.items() if pres}

    ### COHERENCE ###
    # Process coherence
    has_part_dict = coherence.parse_has_part(has_part_file)
    process_coherence, has_part_protein_details = coherence.check_has_part(protein_go_terms, has_part_dict)

    # pathway coherence
    if subontology is None or subontology == 'bp' or subontology == 'mf':
        ec2go_mapping      = coherence.parse_ec2go(ec2go_file)
        pathway_to_go      = coherence.map_pathways_to_go_terms(pathway_file, ec2go_mapping)
        _, metacyc_completed, metacyc_annotated, pathway_details = coherence.analyze_genome(protein_go_terms_ancestors, pathway_to_go, ec2go_mapping)
        if len(metacyc_annotated) > 0:
            metacyc_pct       = (len(metacyc_completed) / len(metacyc_annotated)) * 100
            total_completed   = len(metacyc_completed)
            total_annotated   = len(metacyc_annotated)
            total_incomplete  = total_annotated - total_completed
        else:
            metacyc_pct = 0
            total_completed = 0
            total_annotated = 0
            total_incomplete = 0

    # Protein complex coherence
    if subontology is None or subontology == 'cc':
        # For protein complex coherence, follow the paper: consider only is_a hierarchy
        # (exclude part_of and other relationship edges from ancestor propagation).
        go_is_a = Ontology(ontology_file, with_rels=False)

        protein_go_terms_ancestors_is_a, _ = get_ancestors_with_ontology(protein_go_terms, go_is_a)
        term_to_children, _, _ = coherence.parse_go_ontology(ontology_file, include_part_of=False)
        complex_child_terms = coherence.get_all_child_terms(MACROMOLECULAR_COMPLEX, term_to_children)
        complex_child_terms.add(MACROMOLECULAR_COMPLEX)
        complex_classifications, _ = coherence.classify_complexes(protein_go_terms_ancestors_is_a, complex_child_terms, homodimer_terms=homodimer_terms)
        coherent_count, incoherent_count = coherence.count_complexes(complex_classifications)
        complex_coherence = (coherent_count / (coherent_count + incoherent_count)) * 100 if (coherent_count + incoherent_count) > 0 else 0
        term_names = {t: go.get_term(t)['name'] for t in complex_classifications}

    ### CONSISTENCY ###
    # Taxonomic consistency
    consistency_satisfiable =  check_consistency_fast(protein_go_terms_ancestors, taxa_constraints_file)

    ### OVERVIEW ###
    context = {
        'assembly_name': assembly_name,
        'complete_has_part_percentage': round(process_coherence, 2),
        'has_part_data': has_part_protein_details,
        'satisfiable': consistency_satisfiable,
    }
    if subontology is None or subontology == 'bp':
        # essential terms
        context['essential_percentage'] = round(essential_percentage, 2)
        context['found_terms'] = found_terms
        context['go_categories_core'] = go_core
        context['go_categories_periph'] = go_periph
    if subontology is None or subontology == 'bp' or subontology == 'mf':
        # pathway coherence
        context['ec2go_mapping'] = ec2go_mapping
        context['metacyc_complete_percentage'] = round(metacyc_pct, 2)
        context['metacyc_completed'] = total_completed
        context['metacyc_annotated'] = total_annotated
        context['metacyc_incomplete'] = total_incomplete
        context['pathway_details'] = pathway_details
    if subontology is None or subontology == 'cc':
        # protein complex coherence
        context['complex_coherence'] = round(complex_coherence, 2)
        context['complex_classifications'] = complex_classifications
        context['term_names'] = term_names

    ### INFORMATION CONTENT ###
    ic_depth, ic_breadth, normalized_ic_breadth = compute_ic_depth_breadth_fast(protein_go_terms_specific, go)
    context["ic_depth"] = ic_depth
    context["ic_breadth"] = ic_breadth
    context["normalized_ic_breadth"] = normalized_ic_breadth
    ### SAVE JSON ###
    json_context = context.copy()
    json_context.pop("ec2go_mapping", None)
    json_context.pop("term_names", None)

    if output_file is None:
        output_file = assembly_name + "_report.json"

    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(json_context, f, indent=2, default=str)

    return json_context


def GAEF_evaluation(assembly_name, predictions_file, GAEF_dir, go_file, groovy_flag=False, python_consistency=False, python_ic=True, threshold=0.5, subontology=None, constraints_dir=None):
    """
    Compute GAEF metrics from a predictions file. Thin wrapper around GAEF_evaluation_from_terms().
    """
    protein_go_terms = load_predictions(predictions_file, threshold)
    return GAEF_evaluation_from_terms(assembly_name, protein_go_terms, GAEF_dir, go_file, subontology, constraints_dir=constraints_dir)



def propagate_predictions_dict(predictions_dict, go, ont):
    """
    Propagate prediction scores to ancestors for each protein.
    For each term with score, sets score for all ancestors to max(current, score).

    Args:
        predictions_dict: protein_id -> {go_term: score}
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')

    Returns:
        protein_id -> {go_term: score} with ancestors propagated
    """
    namespace_set = NAMESPACES[ont]
    ancestor_cache = {}
    result = defaultdict(dict)
    for prot_id, term_scores in predictions_dict.items():
        for go_id, score in term_scores.items():
            if not go.has_term(go_id) or go.get_namespace(go_id) != namespace_set:
                continue
            if go_id not in ancestor_cache:
                ancestor_cache[go_id] = go.get_ancestors(go_id)
            result[prot_id][go_id] = max(result[prot_id].get(go_id, 0), score)
            for anc in ancestor_cache[go_id]:
                result[prot_id][anc] = max(result[prot_id].get(anc, 0), score)
    return dict(result)


def load_predictions_prop(predictions_file, go, ont, annotated_proteins=None):
    """
    Load predictions from TSV file and propagate scores to ancestors.
    Optimized with ancestor caching and batch processing.
    
    Args:
        predictions_file: Path to TSV file with format: protein_id\tgo_term1|score1\tgo_term2|score2\t...
                         Lines with only protein_id (no predictions) are also included.
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')
        annotated_proteins: Optional set of protein IDs to filter predictions.
                           Only predictions for proteins in this set will be loaded.
                           This reduces memory usage when evaluating against a specific annotation set.
    
    Returns:
        Dictionary mapping protein_id -> {go_term: score}
        Proteins with no predictions will have an empty dict {}
    """
    predictions = defaultdict(dict)
    # Cache ancestors to avoid recomputation
    ancestor_cache = {}
    namespace_set = NAMESPACES[ont]

    
    with open(predictions_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            it = line.split('\t')
            if len(it) < 1:
                continue
            prot_id = it[0]
            
            # Skip proteins not in annotation set (optimization)
            if annotated_proteins is not None and prot_id not in annotated_proteins:
                continue
            
            preds = it[1:] if len(it) > 1 else []
            term_scores = {}
            # Process all predictions for this protein (may be empty for proteins with no predictions)
            for pred in preds:
                # Parse once
                parts = pred.split('|', 1)
                if len(parts) != 2:
                    continue
                go_id, score_str = parts
                score = float(score_str)

                # Check namespace and term validity
                if not go.has_term(go_id) or go.get_namespace(go_id) != namespace_set:
                    continue
                
                term_scores[go_id] = score
            
            # Batch propagate scores to ancestors (only if there are predictions)
            for go_id, score in term_scores.items():
                # Get ancestors (use cache)
                if go_id not in ancestor_cache:
                    ancestor_cache[go_id] = go.get_ancestors(go_id)
                ancestors = ancestor_cache[go_id]
                
                # Update scores for term and all ancestors
                predictions[prot_id][go_id] = max(predictions[prot_id].get(go_id, 0), score)
                for anc in ancestors:
                    predictions[prot_id][anc] = max(predictions[prot_id].get(anc, 0), score)
    
    return predictions


def load_annotations(annotations_file):
    """
    Load annotations from TSV file. Assumes annotations are propagated to ancestors.
    
    Args:
        annotations_file: Path to TSV file with format: protein_id\tgo_term1\tgo_term2\t...
    
    Returns:
        Dictionary mapping protein_id -> set(go_terms)
    """
    annotations = {}
    if annotations_file is None:
        return annotations
    with open(annotations_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            it = line.split('\t')
            if len(it) < 2:
                continue
            prot_id = it[0]  # expect <prot_id>_<TAXON>
            go_terms = it[1:]
            annotations[prot_id] = set(go_terms)
    return annotations


def evaluate_prediction_metrics(go, ont, annotations, predictions):
    """
    Evaluate prediction metrics for a single ontology namespace following CAFA-style evaluation.

    IMPORTANT: This function follows CAFA evaluation standards where ALL proteins in the
    ground truth annotations are evaluated, including those with no predictions. Proteins
    without predictions are assigned zero scores for all GO terms, which correctly penalizes
    recall while not affecting precision.

    OPTIMIZED VERSION:
    - Pre-allocates eval_preds and labels_matrix in single pass
    - Uses compute_metrics_optimized with IC caching and set intersection
    - Removes redundant precision/recall computation (returned by optimized function)
    - Removes AUC pre-check (handled gracefully in optimized function)

    Args:
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')
        annotations: Dictionary mapping protein_id -> set(go_terms) (with ancestors)
        predictions: Dictionary mapping protein_id -> {go_term: score} (with ancestors)

    Returns:
        Dictionary with metrics: fmax, smin, tmax, wfmax, wtmax, avg_auc, aupr, avgic, precision, recall
        Returns None if no proteins with annotations are found or other errors occur.
    """
    # Pre-compute namespace-filtered GO terms once
    go_set_filtered = go.get_namespace_terms(NAMESPACES[ont])
    root_term = FUNC_DICT[ont]
    if root_term in go_set_filtered:
        go_set_filtered.remove(root_term)

    # Get all terms from annotations, filtered by namespace (use set intersection)
    terms = set()
    for go_terms in annotations.values():
        terms.update(go_terms & go_set_filtered)
    terms.discard(root_term)

    if len(terms) == 0:
        return None

    terms = list(terms)
    terms_dict = {v: i for i, v in enumerate(terms)}
    go_set = set(terms_dict)
    num_terms = len(terms)

    # First pass: count valid proteins to pre-allocate arrays
    ancestor_cache = {}
    valid_proteins = []

    for prot_id, annot_terms in annotations.items():
        annots = go_set & annot_terms  # Set intersection
        if len(annots) > 0:
            valid_proteins.append((prot_id, annot_terms))

    if len(valid_proteins) == 0:
        return None

    num_proteins = len(valid_proteins)

    # Pre-allocate matrices (single allocation instead of list of arrays)
    eval_preds = np.zeros((num_proteins, num_terms), dtype=np.float32)
    labels_matrix = np.zeros((num_proteins, num_terms), dtype=np.float32)

    # Single pass: build proteins list, prop_annotations, and fill matrices
    proteins = []
    prop_annotations = []

    for idx, (prot_id, annot_terms) in enumerate(valid_proteins):
        annots = go_set & annot_terms

        # Get propagated annotations (all ancestors) - use cache
        prop_annots = set()
        for go_id in annots:
            if go_id not in ancestor_cache:
                ancestor_cache[go_id] = go.get_ancestors(go_id)
            prop_annots |= ancestor_cache[go_id]
        # Filter propagated annotations by namespace (use set intersection)
        prop_annots = prop_annots & go_set_filtered

        proteins.append(prot_id)
        prop_annotations.append(prop_annots)

        # Fill labels_matrix row
        for go_id in prop_annots:
            if go_id in terms_dict:
                labels_matrix[idx, terms_dict[go_id]] = 1

        # Fill eval_preds row directly (no intermediate array creation)
        if prot_id in predictions:
            scores = predictions[prot_id]
            for go_id, score in scores.items():
                if go_id in terms_dict:
                    eval_preds[idx, terms_dict[go_id]] = score

    # Create test dataframe
    test_df = pd.DataFrame({
        'proteins': proteins,
        'prop_annotations': prop_annotations,
        'exp_annotations': prop_annotations
    })

    # Calculate IC for terms
    try:
        go.calculate_ic(prop_annotations)
    except Exception as e:
        print(f"[Error] Error calculating IC for {ont}: {e}")
        return None

    # CAFA-style: Check if there are NO predictions at all
    if np.max(eval_preds) == 0:
        print(f"[Warning] No predictions found for any annotated protein in {ont}. Assigning zero metrics (CAFA-style).")
        return {
            'fmax': 0.0,
            'smin': float('inf'),
            'tmax': 0.0,
            'wfmax': 0.0,
            'wtmax': 0.0,
            'avg_auc': 0.0,
            'aupr': 0.0,
            'avgic': 0.0,
            'precision': 0.0,
            'recall': 0.0,
            'sensitivity': 0.0,
            'specificity': 0.0,
            'fmax_spec_match': 0.0,
            'tp': 0,
            'fp': 0,
            'fn': 0
        }

    # Compute metrics using optimized function (passes pre-built labels_matrix)
    try:
        result = compute_metrics_optimized(
            test_df, go, terms_dict, terms, ont, eval_preds, labels_matrix=labels_matrix)

        fmax, smin, tmax, wfmax, wtmax, avg_auc, aupr, avgic, fmax_spec_match, precision, recall, tp, fp, fn, sensitivity, specificity = result

        return {
            'fmax': fmax,
            'smin': smin,
            'tmax': tmax,
            'wfmax': wfmax,
            'wtmax': wtmax,
            'avg_auc': avg_auc,
            'aupr': aupr,
            'avgic': avgic,
            'precision': precision,
            'recall': recall,
            'sensitivity': sensitivity,
            'specificity': specificity,
            'fmax_spec_match': fmax_spec_match,
            'tp': tp,
            'fp': fp,
            'fn': fn
        }
    except ZeroDivisionError as e:
        print(f"[Error] Division by zero in AUC computation for {ont}: {e}")
        print(f"  This happens when no GO terms have valid distributions for per-term AUC.")
        print(f"  Computing F-max and other metrics without AUC...")

        # Fallback: compute metrics with AUC set to 0
        # Use optimized evaluate_annotations with IC cache
        go_set_for_labels = go.get_namespace_terms(NAMESPACES[ont])
        if root_term in go_set_for_labels:
            go_set_for_labels.remove(root_term)

        # Use set intersection instead of map/lambda/filter
        labels = [annots & go_set_for_labels for annots in prop_annotations]
        spec_labels = labels  # Same as prop_annotations for our case

        fmax = 0.0
        tmax = 0.0
        wfmax = 0.0
        wtmax = 0.0
        avgic = 0.0
        smin = 1000000.0
        fmax_spec_match = 0
        precisions = []
        recalls = []
        precision_at_tmax = 0.0
        recall_at_tmax = 0.0
        tp_at_tmax = 0
        fp_at_tmax = 0
        fn_at_tmax = 0

        # IC cache for optimized evaluation
        ic_cache = {}

        # Threshold sweep
        for t in range(0, 101):
            threshold = t / 100.0

            # Vectorized threshold + set intersection
            above_threshold = eval_preds >= threshold
            preds_at_t = []
            for i in range(num_proteins):
                above_indices = np.where(above_threshold[i])[0]
                annots = {terms[j] for j in above_indices} & go_set_for_labels
                preds_at_t.append(annots)

            # Evaluate at this threshold
            fscore, prec, rec, s, ru, mi, fps, fns, avg_ic, wf = evaluate_annotations_optimized(
                go, labels, preds_at_t, ic_cache)

            spec_match = sum(len(spec_labels[i] & preds_at_t[i]) for i in range(num_proteins))

            tp_sum = sum(len(labels[i] & preds_at_t[i]) for i in range(num_proteins) if len(labels[i]) > 0)
            fp_sum = sum(len(fp) for fp in fps)
            fn_sum = sum(len(fn) for fn in fns)

            precisions.append(prec)
            recalls.append(rec)

            if fmax < fscore:
                fmax = fscore
                tmax = threshold
                avgic = avg_ic
                fmax_spec_match = spec_match
                precision_at_tmax = prec
                recall_at_tmax = rec
                tp_at_tmax = tp_sum
                fp_at_tmax = fp_sum
                fn_at_tmax = fn_sum
            if wfmax < wf:
                wfmax = wf
                wtmax = threshold
            if smin > s:
                smin = s

        # Compute AUPR
        precisions = np.array(precisions)
        recalls = np.array(recalls)
        sorted_index = np.argsort(recalls)
        recalls = recalls[sorted_index]
        precisions = precisions[sorted_index]
        aupr = np.trapz(precisions, recalls)

        # Term-centric sensitivity and specificity at F-max threshold
        above_at_tmax = eval_preds >= tmax
        sens_list, spec_list = [], []
        for i in range(num_terms):
            labels_i = labels_matrix[:, i]
            preds_i = above_at_tmax[:, i]
            tp_t = np.sum((labels_i == 1) & preds_i)
            fn_t = np.sum((labels_i == 1) & ~preds_i)
            fp_t = np.sum((labels_i == 0) & preds_i)
            tn_t = np.sum((labels_i == 0) & ~preds_i)
            pos_n = tp_t + fn_t
            neg_n = tn_t + fp_t
            if pos_n > 0:
                sens_list.append(tp_t / pos_n)
            if neg_n > 0:
                spec_list.append(tn_t / neg_n)
        sensitivity_at_tmax = float(np.mean(sens_list)) if sens_list else 0.0
        specificity_at_tmax = float(np.mean(spec_list)) if spec_list else 0.0

        avg_auc = 0.0
        print(f"  Successfully computed F-max={fmax:.4f} without AUC")

        return {
            'fmax': fmax,
            'smin': smin,
            'tmax': tmax,
            'wfmax': wfmax,
            'wtmax': wtmax,
            'avg_auc': avg_auc,
            'aupr': aupr,
            'avgic': avgic,
            'precision': precision_at_tmax,
            'recall': recall_at_tmax,
            'sensitivity': sensitivity_at_tmax,
            'specificity': specificity_at_tmax,
            'fmax_spec_match': fmax_spec_match,
            'tp': tp_at_tmax,
            'fp': fp_at_tmax,
            'fn': fn_at_tmax
        }
    except Exception as e:
        print(f"[Error] Error computing metrics for {ont}: {e}")
        import traceback
        traceback.print_exc()
        return None


def _build_term_matrices(go, ont, annotations, predictions):
    """Build aligned eval_preds and labels_matrix for per-term analysis.

    Args:
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')
        annotations: Dictionary mapping protein_id -> set(go_terms) (with ancestors)
        predictions: Dictionary mapping protein_id -> {go_term: score} (with ancestors, or binary 0/1)

    Returns:
        Tuple of (terms, terms_dict, eval_preds, labels_matrix, valid_proteins, prop_annotations)
        or (None, None, None, None, None, None) on error.
    """
    go_set_filtered = go.get_namespace_terms(NAMESPACES[ont])
    root_term = FUNC_DICT[ont]
    if root_term in go_set_filtered:
        go_set_filtered.remove(root_term)

    terms = set()
    for go_terms in annotations.values():
        terms.update(go_terms & go_set_filtered)
    terms.discard(root_term)

    if len(terms) == 0:
        return None, None, None, None, None, None

    terms = list(terms)
    terms_dict = {v: i for i, v in enumerate(terms)}
    go_set = set(terms_dict)

    valid_proteins = []
    for prot_id, annot_terms in annotations.items():
        annots = go_set & annot_terms
        if len(annots) > 0:
            valid_proteins.append((prot_id, annot_terms))

    if len(valid_proteins) == 0:
        return None, None, None, None, None, None

    num_proteins = len(valid_proteins)
    num_terms = len(terms)
    eval_preds = np.zeros((num_proteins, num_terms), dtype=np.float32)
    labels_matrix = np.zeros((num_proteins, num_terms), dtype=np.float32)
    prop_annotations = []
    ancestor_cache = {}

    for idx, (prot_id, annot_terms) in enumerate(valid_proteins):
        annots = go_set & annot_terms
        prop_annots = set()
        for go_id in annots:
            if go_id not in ancestor_cache:
                ancestor_cache[go_id] = go.get_ancestors(go_id)
            prop_annots |= ancestor_cache[go_id]
        prop_annots = prop_annots & go_set_filtered
        prop_annotations.append(prop_annots)

        for go_id in prop_annots:
            if go_id in terms_dict:
                labels_matrix[idx, terms_dict[go_id]] = 1

        if prot_id in predictions:
            scores = predictions[prot_id]
            for go_id, score in scores.items():
                if go_id in terms_dict:
                    eval_preds[idx, terms_dict[go_id]] = score

    try:
        go.calculate_ic(prop_annotations)
    except Exception as e:
        print(f"[Error] Error calculating IC for {ont}: {e}")
        return None, None, None, None, None, None

    return terms, terms_dict, eval_preds, labels_matrix, valid_proteins, prop_annotations


def compute_per_term_metrics(go, ont, annotations, predictions, threshold=None):
    """
    Compute per-term sensitivity, specificity, and n_annotations.

    Used for comparing term-specific metrics between predictions and optimized annotations.
    For scored predictions, threshold=None finds tmax via F-max. For binary (e.g. optimized),
    pass threshold=0.5 to treat any score > 0.5 as predicted.

    Args:
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')
        annotations: Dictionary mapping protein_id -> set(go_terms) (with ancestors)
        predictions: Dictionary mapping protein_id -> {go_term: score} (with ancestors, or binary 0/1)
        threshold: Optional. If None, runs threshold sweep to find tmax. If float (e.g. 0.5), uses that.

    Returns:
        Tuple of (per_term_dict, tmax) or (None, None) on error.
        per_term_dict: {term_id: {'sensitivity': float, 'specificity': float, 'n_annotations': int}}
    """
    result = _build_term_matrices(go, ont, annotations, predictions)
    terms, terms_dict, eval_preds, labels_matrix, valid_proteins, prop_annotations = result
    if terms is None:
        return None, None

    if threshold is None:
        test_df = pd.DataFrame({
            'proteins': [p[0] for p in valid_proteins],
            'prop_annotations': prop_annotations,
            'exp_annotations': prop_annotations
        })
        try:
            result = compute_metrics_optimized(
                test_df, go, terms_dict, terms, ont, eval_preds, labels_matrix=labels_matrix)
            tmax = result[2]
        except ZeroDivisionError:
            tmax = 0.5
    else:
        tmax = threshold

    above_threshold = eval_preds >= tmax
    per_term = {}
    for i, term_id in enumerate(terms):
        labels_i = labels_matrix[:, i]
        preds_i = above_threshold[:, i]
        tp = int(np.sum((labels_i == 1) & preds_i))
        fn = int(np.sum((labels_i == 1) & ~preds_i))
        fp = int(np.sum((labels_i == 0) & preds_i))
        tn = int(np.sum((labels_i == 0) & ~preds_i))
        n_annotations = tp + fn
        pos_n = tp + fn
        neg_n = tn + fp
        sens = tp / pos_n if pos_n > 0 else 0.0
        spec = tn / neg_n if neg_n > 0 else 0.0
        per_term[term_id] = {
            'sensitivity': float(sens),
            'specificity': float(spec),
            'n_annotations': n_annotations
        }

    return per_term, tmax


def compute_per_term_auc(go, ont, annotations, predictions):
    """Compute per-term ROC AUC from continuous prediction scores.

    Args:
        go: Ontology object
        ont: Ontology namespace ('cc', 'bp', or 'mf')
        annotations: Dictionary mapping protein_id -> set(go_terms) (with ancestors)
        predictions: Dictionary mapping protein_id -> {go_term: score} (continuous scores)

    Returns:
        Dict of {term_id: {'auc': float, 'n_annotations': int}} or None on error.
        AUC is NaN for terms with no positive or no negative examples.
    """
    result = _build_term_matrices(go, ont, annotations, predictions)
    terms, terms_dict, eval_preds, labels_matrix, valid_proteins, prop_annotations = result
    if terms is None:
        return None

    num_proteins = len(valid_proteins)
    per_term_auc = {}
    for i, term_id in enumerate(terms):
        pos_n = int(np.sum(labels_matrix[:, i]))
        neg_n = num_proteins - pos_n
        n_annotations = pos_n
        if pos_n > 0 and neg_n > 0:
            auc_val = float(compute_roc_auc(labels_matrix[:, i], eval_preds[:, i]))
        else:
            auc_val = float('nan')
        per_term_auc[term_id] = {
            'auc': auc_val,
            'n_annotations': n_annotations,
        }

    return per_term_auc


def main(assembly_name, annotations_file, predictions_file, GAEF_dir, go_file, groovy_flag=False, python_consistency=False, python_ic=True, threshold=0.5, output_file=None, subontology=None):
    # Run GAEF evaluation
    json_context = GAEF_evaluation(assembly_name, predictions_file, GAEF_dir, go_file, groovy_flag, python_consistency, python_ic, threshold, subontology)
    
    print(f"--- GAEF evaluation for {assembly_name} (threshold: {threshold}) ---")
    if subontology is None or subontology == 'bp':
        print(f"- Essential percentage: {json_context['essential_percentage']}")
    print(f"- Taxonomic consistency: {json_context['satisfiable']}")
    print(f"- Process coherence: {json_context['complete_has_part_percentage']}")
    if subontology is None or subontology == 'cc':
        print(f"- Complex coherence: {json_context['complex_coherence']}")
    if subontology is None or subontology == 'bp' or subontology == 'mf':
        print(f"- Metacyc complete percentage: {json_context['metacyc_complete_percentage']}")
    
    # Evaluate prediction metrics for each ontology namespace
    ontology_file = go_file
    go = Ontology(ontology_file)
    # Load annotations once and reuse across all namespaces
    all_annotations = load_annotations(annotations_file)
    
    prediction_metrics = {}
    try:
        # Determine which ontologies to evaluate
        if subontology is not None:
            ontologies_to_evaluate = [subontology]
        else:
            ontologies_to_evaluate = ['cc', 'bp', 'mf']
        
        # Evaluate metrics for each ontology namespace
        for ont in ontologies_to_evaluate:
            print(f"\nProcessing {ont} ontology for {assembly_name} (threshold: {threshold})...")
            
            # Reuse loaded annotations (already filtered by namespace in evaluate_prediction_metrics)
            if len(all_annotations) == 0:
                print(f"[Warning] No annotations found for {ont} namespace")
                continue
            
            # Load predictions for this namespace (only for annotated proteins - optimization)
            annotated_protein_ids = set(all_annotations.keys())
            predictions = load_predictions_prop(predictions_file, go, ont, annotated_proteins=annotated_protein_ids)
            
            print(f"Loaded predictions for {len(predictions)} proteins (where {len(annotated_protein_ids)} are annotated)")
            
            # Note: Even if len(predictions) == 0, we still evaluate (CAFA-style)
            # Proteins without predictions will receive zero scores for all terms
            # This correctly penalizes models that fail to make predictions
            
            # Evaluate metrics
            metrics = evaluate_prediction_metrics(go, ont, all_annotations, predictions)
            if metrics is not None:
                prediction_metrics[ont] = metrics
                print(f"---- {ont.upper()} ----")
                print(f"Fmax: {metrics['fmax']:.5f}, Smin: {metrics['smin']:.5f}, Threshold: {metrics['tmax']:.5f}, WFmax: {metrics['wfmax']:.5f}, WTmax: {metrics['wtmax']:.5f}")
                print(f"AUC: {metrics['avg_auc']:.5f}")
                print(f"AUPR: {metrics['aupr']:.5f}")
                print(f"AVGIC: {metrics['avgic']:.5f}")
                print(f"Precision: {metrics['precision']:.5f}")
                print(f"Recall: {metrics['recall']:.5f}")
                print(f"Sensitivity (term-centric): {metrics.get('sensitivity', 0.0):.5f}")
                print(f"Specificity (term-centric): {metrics.get('specificity', 0.0):.5f}")
                print(f"TP: {metrics.get('tp', 'N/A')}, FP: {metrics.get('fp', 'N/A')}, FN: {metrics.get('fn', 'N/A')}")
            else:
                print(f"[Warning] Could not compute metrics for {ont} namespace (possible IC calculation failure or insufficient data)")
        
        # Calculate overall metrics (average across namespaces)
        # If only one subontology is evaluated, overall metrics are the same as that subontology
        if len(prediction_metrics) > 0:
            overall_metrics = {}
            metric_keys = ['fmax', 'smin', 'wfmax', 'avg_auc', 'aupr', 'avgic', 'precision', 'recall', 'sensitivity', 'specificity', 'tp', 'fp', 'fn']
            for key in metric_keys:
                values = [prediction_metrics[ont][key] for ont in prediction_metrics.keys() if key in prediction_metrics[ont]]
                if values:
                    overall_metrics[f'avg_{key}'] = sum(values) / len(values)
            
            prediction_metrics['overall'] = overall_metrics
            json_context['prediction_metrics'] = prediction_metrics
            
            # Update the JSON file with prediction metrics
            json_context_copy = json_context.copy()

            if output_file is None:
                output_file = f"{assembly_name}_report.json"

            # create the directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            with open(output_file, "w") as f:
                json.dump(json_context_copy, f, indent=2, default=str)
            
            print(f"\nPrediction metrics computed for {len(prediction_metrics) - 1} ontology namespace(s)")
            print(f"Updated {output_file} with prediction metrics")
        else:
            print("[Warning] No prediction metrics were computed")
    except Exception as e:
        print(f"[Error] Error evaluating prediction metrics: {e}")
        import traceback
        traceback.print_exc()
    
    return json_context



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate GAEF metrics')
    parser.add_argument('--assembly_name', required=True, help='Name of the assembly/genome')
    parser.add_argument('--annotations_file', required=False, default=None, help='Path to GAEF format annotation file')
    parser.add_argument('--predictions_file', required=True, help='Path to predictions file (TSV format: protein_id\\tgo_term1|score1\\tgo_term2|score2\\t...). If provided, prediction metrics will be computed.')
    parser.add_argument('--GAEF_dir', default="../../GAEF", help='Path to GAEF directory')
    parser.add_argument('--go_file', default="data/go-basic.obo", help='Path to GO file')  # 2025-10 version
    parser.add_argument('--threshold', type=float, default=0.5, help='Threshold for GO term scores')
    parser.add_argument('--groovy_flag', action='store_true', help='Use Groovy scripts for IC calculation and taxonomic consistency')
    parser.add_argument('--python_consistency', action='store_true', help='Use Python implementation for taxonomic consistency (faster, no JVM)')
    parser.add_argument('--python_ic', action='store_true', help='Use Python implementation for IC calculation (default: True)')
    parser.add_argument('--output_file', default=None, help='Path to the output file')
    parser.add_argument('--subontology', type=str, default=None, choices=['cc', 'bp', 'mf'],
                       help='Subontology to evaluate (default: None, evaluates all)')
    args = parser.parse_args()
    main(
        assembly_name=args.assembly_name, 
        annotations_file=args.annotations_file, 
        predictions_file=args.predictions_file, 
        GAEF_dir=args.GAEF_dir, 
        go_file=args.go_file,
        groovy_flag=args.groovy_flag, 
        python_consistency=args.python_consistency, 
        python_ic=args.python_ic, 
        threshold=args.threshold,
        output_file=args.output_file,
        subontology=args.subontology
        )