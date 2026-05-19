#!/usr/bin/env python3
"""
Heuristic ablation study for the complex coherence CP-SAT solver.

Evaluates eight heuristic configurations on a single organism's CC predictions:

  Config       H1(incoherent)  H2(sparse)  H3(top-k)
  ───────────  ──────────────  ──────────  ─────────
  Naive              -              -           -
  +H1                ✓              -           -
  +H2                -              ✓           -
  +H3                -              -           ✓
  +H1+H3             ✓              -           ✓
  +H2+H3             -              ✓           ✓
  +H1+H2             ✓              ✓           -
  +H1+H2+H3          ✓              ✓           ✓

With --dry_run, only variable/constraint counts are computed (no CP-SAT model
is built), so OOM/timeout cannot occur.

Metrics reported per config:
  n_vars        number of Boolean variables (would be) created
  n_constraints number of constraints (would be) added
  setup_s       model construction time in seconds (0 in dry-run)
  solve_s       solver wall-clock time in seconds (0 in dry-run)
  flips         number of annotation changes in the solution (- in dry-run)
  status        OPTIMAL | TIMEOUT | OOM | INFEASIBLE | DRY_RUN

Usage:
  python heuristic_ablation.py \\
      --predictions <path> --complexes <path> --go_hierarchy <path> \\
      --taxon_id 272844 --threshold 0.3 --top_k 50 --timeout 300

  python heuristic_ablation.py --dry_run \\
      --predictions <path> --complexes <path> --go_hierarchy <path> \\
      --taxon_id 272844 --threshold 0.3 --top_k 50
"""

import argparse
import csv
import sys
import time
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from ortools.sat.python import cp_model

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ROOT_TERMS = {"GO:0005575", "GO:0008150", "GO:0003674"}
MACROMOLECULAR_COMPLEX = "GO:0032991"


# ---------------------------------------------------------------------------
# Data loading  (unchanged from adjust_ortools.py)
# ---------------------------------------------------------------------------

def load_predictions(file_path: str) -> Dict[str, Dict[str, float]]:
    predictions = defaultdict(dict)
    with open(file_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) < 2:
                continue
            protein = row[0]
            predictions[protein] = {
                go: float(score)
                for field in row[1:]
                for go, score in [field.split("|")]
            }
    return predictions


def load_homodimer_terms(file_path: str) -> Set[str]:
    terms = set()
    with open(file_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[0] == "GO_term":
                continue
            if row[1] in ("h", "a"):
                terms.add(row[0])
    return terms


def load_go_hierarchy(file_path: str) -> Dict[str, Set[str]]:
    hierarchy = defaultdict(set)
    with open(file_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            hierarchy[row[0]].add(row[1])   # child -> parents
    return hierarchy


# ---------------------------------------------------------------------------
# Derived data structures
# ---------------------------------------------------------------------------

def get_all_children(go_term: str, hierarchy: Dict[str, Set[str]]) -> Set[str]:
    children, queue = set(), [go_term]
    while queue:
        cur = queue.pop(0)
        for child in sorted(c for c in hierarchy if cur in hierarchy[c]):
            if child not in children:
                children.add(child)
                queue.append(child)
    return children


def get_heteromeric_complexes(homodimer_terms: Set[str],
                              hierarchy: Dict[str, Set[str]]) -> Set[str]:
    all_complex = get_all_children(MACROMOLECULAR_COMPLEX, hierarchy)
    all_complex.add(MACROMOLECULAR_COMPLEX)
    return all_complex - homodimer_terms


def get_annotated_predictions(predictions: Dict[str, Dict[str, float]],
                               threshold: float) -> Tuple[Dict[str, Set[str]], Set[str]]:
    annotated = defaultdict(set, {
        p: {g for g, s in scores.items() if s > threshold}
        for p, scores in predictions.items()
    })
    all_terms = set().union(*annotated.values())
    return annotated, all_terms


def get_incoherent_complexes(annotated: Dict[str, Set[str]],
                              complex_terms: Set[str]) -> Set[str]:
    """Return the subset of complex_terms predicted for exactly one protein."""
    incoherent = set()
    for term in complex_terms:
        count = sum(1 for p in annotated if term in annotated[p])
        if count == 1:
            incoherent.add(term)
    return incoherent


def compute_flip_costs(predictions: Dict[str, Dict[str, float]],
                       threshold: float) -> Dict[str, Dict[str, float]]:
    costs: Dict[str, Dict[str, float]] = defaultdict(
        lambda: defaultdict(lambda: threshold)
    )
    for protein, scores in predictions.items():
        for go, score in scores.items():
            costs[protein][go] = abs(score - threshold)
    return costs


def get_ancestors(go_term: str, hierarchy: Dict[str, Set[str]]) -> Set[str]:
    if go_term in ROOT_TERMS:
        return set()
    ancestors = set()
    for parent in hierarchy[go_term]:
        ancestors.add(parent)
        ancestors.update(get_ancestors(parent, hierarchy))
    return ancestors


def get_descendants(go_term: str, hierarchy: Dict[str, Set[str]]) -> Set[str]:
    desc = set()
    for child in [c for c in hierarchy if go_term in hierarchy[c]]:
        desc.add(child)
        desc.update(get_descendants(child, hierarchy))
    return desc


def get_participating_proteins(annotated: Dict[str, Set[str]],
                                active_terms: Set[str],
                                predictions: Dict[str, Dict[str, float]],
                                top_k: int = None) -> Dict[str, List[str]]:
    """
    For each active complex term, return the list of candidate proteins.
    With top_k: ranked by prediction score, truncated to top_k.
    Without top_k: all proteins in the proteome.
    """
    result: Dict[str, List[str]] = {}
    for term in active_terms:
        scored = [(p, predictions[p].get(term, 0.0)) for p in predictions]
        scored.sort(key=lambda x: x[1], reverse=True)
        result[term] = [p for p, _ in (scored[:top_k] if top_k else scored)]
    return result


def get_go_protein_pairs(annotated: Dict[str, Set[str]],
                          hierarchy: Dict[str, Set[str]],
                          participating: Dict[str, List[str]]) -> Set[Tuple[str, str]]:
    """Sparse (H2) variable set: demotion paths for predicted terms,
    promotion paths (ancestors) for non-predicted terms."""
    pairs: Set[Tuple[str, str]] = set()
    for term, proteins in participating.items():
        for protein in proteins:
            pairs.add((term, protein))
            if term in annotated[protein]:
                for desc in get_descendants(term, hierarchy):
                    if desc in annotated[protein]:
                        pairs.add((desc, protein))
            else:
                for anc in get_ancestors(term, hierarchy):
                    if anc not in annotated[protein]:
                        pairs.add((anc, protein))
    return pairs


def get_hierarchy_terms(hierarchy: Dict[str, Set[str]]) -> Set[str]:
    """Return the set of all GO terms mentioned in the hierarchy."""
    terms: Set[str] = set()
    for child, parents in hierarchy.items():
        terms.add(child)
        terms.update(parents)
    return terms


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

class RunResult:
    def __init__(self, config: str):
        self.config = config
        self.n_vars = 0
        self.n_constraints = 0
        self.setup_s = 0.0
        self.solve_s = 0.0
        self.flips = 0
        self.status = "NOT_RUN"


# ---------------------------------------------------------------------------
# Dry-run counting functions
# ---------------------------------------------------------------------------

def count_full_vars_config(annotated: Dict[str, Set[str]],
                           hierarchy: Dict[str, Set[str]],
                           active_complexes: Set[str],
                           config_label: str,
                           participating_proteins: Optional[Dict[str, List[str]]] = None
                           ) -> RunResult:
    """Count variables and constraints for a full-variable-set configuration
    without building a CP-SAT model."""
    result = RunResult(config_label)
    t0 = time.time()

    hierarchy_terms = get_hierarchy_terms(hierarchy)
    n_proteins = len(annotated)
    n_pairs = len(hierarchy_terms) * n_proteins
    result.n_vars = 2 * n_pairs   # go_var + flip_var per pair

    # Hierarchy constraints: each hierarchy edge × each protein
    n_hier_edges = sum(len(parents) for parents in hierarchy.values())
    n_hier = n_hier_edges * n_proteins

    # Complex coherence constraints
    active_in_hierarchy = active_complexes & hierarchy_terms
    n_complex = 0
    for term in active_in_hierarchy:
        if participating_proteins is not None:
            n_cand = len(participating_proteins.get(term, []))
        else:
            n_cand = n_proteins
        if n_cand > 1:
            n_complex += n_cand

    # XOR constraints: 1 per pair
    n_xor = n_pairs

    result.n_constraints = n_hier + n_complex + n_xor
    result.setup_s = time.time() - t0
    result.status = "DRY_RUN"
    return result


def count_sparse_config(annotated: Dict[str, Set[str]],
                        hierarchy: Dict[str, Set[str]],
                        go_protein_pairs: Set[Tuple[str, str]],
                        participating: Dict[str, List[str]],
                        config_label: str) -> RunResult:
    """Count variables and constraints for a sparse-variable-set configuration
    without building a CP-SAT model."""
    result = RunResult(config_label)
    t0 = time.time()

    result.n_vars = 2 * len(go_protein_pairs)

    n_con = 0
    # Hierarchy constraints
    for child, protein in go_protein_pairs:
        for parent in hierarchy.get(child, set()):
            if (parent, protein) in go_protein_pairs:
                n_con += 1

    # Complex coherence: 3 constraints per active term (sum==1, sum!=1, and-not)
    for term, proteins in participating.items():
        if any((term, p) in go_protein_pairs for p in proteins):
            n_con += 3

    # XOR constraints
    n_con += len(go_protein_pairs)

    result.n_constraints = n_con
    result.setup_s = time.time() - t0
    result.status = "DRY_RUN"
    return result


# ---------------------------------------------------------------------------
# Full-variable-set solver  (Naive, +H1, +H3, +H1+H3)
# ---------------------------------------------------------------------------

def solve_full_vars(annotated: Dict[str, Set[str]],
                    hierarchy: Dict[str, Set[str]],
                    flip_costs: Dict[str, Dict[str, float]],
                    active_complexes: Set[str],
                    timeout: int,
                    config_label: str,
                    participating_proteins: Optional[Dict[str, List[str]]] = None
                    ) -> RunResult:
    """
    Solver using the full (hierarchy_term × protein) variable set.

    active_complexes: which complex terms get coherence constraints.
    participating_proteins: if provided, limits which proteins get coherence
        constraints per complex term (H3).  If None, all proteins are used.
    """
    result = RunResult(config_label)
    t0 = time.time()

    model = cp_model.CpModel()

    hierarchy_terms = get_hierarchy_terms(hierarchy)

    all_pairs: Set[Tuple[str, str]] = set()
    for term in sorted(hierarchy_terms):
        for protein in sorted(annotated):
            all_pairs.add((term, protein))

    go_vars: Dict[Tuple[str, str], cp_model.IntVar] = {}
    flip_vars: Dict[Tuple[str, str], cp_model.IntVar] = {}
    for pair in sorted(all_pairs):
        go_vars[pair] = model.NewBoolVar(f"go_{pair[0]}_{pair[1]}")
        flip_vars[pair] = model.NewBoolVar(f"flip_{pair[0]}_{pair[1]}")

    result.n_vars = len(go_vars) + len(flip_vars)
    print(f"    vars={result.n_vars:,}")

    # GO hierarchy constraints
    n_con = 0
    for child, protein in sorted(all_pairs):
        for parent in sorted(hierarchy[child]):
            if (parent, protein) in go_vars:
                model.AddImplication(go_vars[(child, protein)],
                                     go_vars[(parent, protein)])
                n_con += 1

    # Complex coherence constraints
    for term in sorted(active_complexes):
        if participating_proteins is not None:
            candidates = participating_proteins.get(term, [])
        else:
            candidates = sorted(annotated)
        for protein in candidates:
            if (term, protein) in go_vars:
                others = [go_vars[(term, p)] for p in candidates
                          if p != protein and (term, p) in go_vars]
                if others:
                    model.AddBoolOr(others).OnlyEnforceIf(go_vars[(term, protein)])
                    n_con += 1

    # XOR constraints (annotated: go+flip=1; non-annotated: go=flip)
    for term, protein in sorted(all_pairs):
        if term in annotated[protein]:
            model.Add(go_vars[(term, protein)] + flip_vars[(term, protein)] == 1)
        else:
            model.Add(go_vars[(term, protein)] == flip_vars[(term, protein)])
        n_con += 1

    # Objective
    cost_scale = 1_000_000
    model.Minimize(sum(
        flip_vars[pair] * max(1, int(flip_costs[pair[1]][pair[0]] * cost_scale))
        for pair in sorted(all_pairs)
    ))

    result.n_constraints = n_con
    result.setup_s = time.time() - t0
    print(f"    constraints={result.n_constraints:,}")
    print(f"    setup={result.setup_s:.2f}s")
    if result.setup_s > timeout:
        result.status = "TIMEOUT"
        print(f"    solver is skipped due to setup timeout")
        return result
    # Solve with timeout
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    t1 = time.time()
    status = solver.Solve(model)
    result.solve_s = time.time() - t1
    print(f"    solve={result.solve_s:.2f}s")
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result.status = "OPTIMAL" if status == cp_model.OPTIMAL else "FEASIBLE"
        result.flips = sum(
            1 for pair in all_pairs
            if pair in flip_vars and solver.BooleanValue(flip_vars[pair])
        )
    elif status == cp_model.UNKNOWN:
        result.status = "TIMEOUT"
    else:
        result.status = "INFEASIBLE"
    print(f"    flips={result.flips}  status={result.status}")
    return result


# ---------------------------------------------------------------------------
# Sparse-variable-set solver  (+H2, +H1+H2, +H1+H2+H3, +H2+H3)
# ---------------------------------------------------------------------------

def solve_sparse_vars(annotated: Dict[str, Set[str]],
                      hierarchy: Dict[str, Set[str]],
                      flip_costs: Dict[str, Dict[str, float]],
                      participating: Dict[str, List[str]],
                      go_protein_pairs: Set[Tuple[str, str]],
                      timeout: int,
                      config_label: str) -> RunResult:
    """
    Solver using the sparse (H2) variable set.
    participating controls candidate proteins per complex (all or top-k).
    """
    result = RunResult(config_label)
    t0 = time.time()

    model = cp_model.CpModel()

    pairs = sorted(go_protein_pairs)
    go_vars: Dict[Tuple[str, str], cp_model.IntVar] = {}
    flip_vars: Dict[Tuple[str, str], cp_model.IntVar] = {}
    for pair in pairs:
        go_vars[pair] = model.NewBoolVar(f"go_{pair[0]}_{pair[1]}")
        flip_vars[pair] = model.NewBoolVar(f"flip_{pair[0]}_{pair[1]}")

    result.n_vars = len(go_vars) + len(flip_vars)
    print(f"    vars={result.n_vars:,}")
    n_con = 0

    # GO hierarchy constraints
    for child, protein in pairs:
        for parent in sorted(hierarchy[child]):
            if (parent, protein) in go_vars:
                model.AddImplication(go_vars[(child, protein)],
                                     go_vars[(parent, protein)])
                n_con += 1

    # Complex coherence constraints
    for term, proteins in participating.items():
        complex_vars = [go_vars[(term, p)] for p in proteins
                        if (term, p) in go_vars]
        if complex_vars:
            sum_is_one = model.NewBoolVar(f"sum_is_one_{term}")
            model.Add(sum(complex_vars) == 1).OnlyEnforceIf(sum_is_one)
            model.Add(sum(complex_vars) != 1).OnlyEnforceIf(sum_is_one.Not())
            model.AddBoolAnd([sum_is_one.Not()])
            n_con += 3

    # XOR constraints
    for term, protein in pairs:
        if term in annotated[protein]:
            model.Add(go_vars[(term, protein)] + flip_vars[(term, protein)] == 1)
        else:
            model.Add(go_vars[(term, protein)] == flip_vars[(term, protein)])
        n_con += 1

    # Objective
    cost_scale = 1_000_000
    model.Minimize(sum(
        flip_vars[pair] * max(1, int(flip_costs[pair[1]][pair[0]] * cost_scale))
        for pair in pairs
    ))

    result.n_constraints = n_con
    result.setup_s = time.time() - t0
    print(f"    constraints={result.n_constraints:,}")
    print(f"    setup={result.setup_s:.2f}s")
    if result.setup_s > timeout:
        result.status = "TIMEOUT"
        print(f"    solver is skipped due to setup timeout")
        return result
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    t1 = time.time()
    status = solver.Solve(model)
    result.solve_s = time.time() - t1
    print(f"    solve={result.solve_s:.2f}s")
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        result.status = "OPTIMAL" if status == cp_model.OPTIMAL else "FEASIBLE"
        result.flips = sum(
            1 for pair in pairs if solver.BooleanValue(flip_vars[pair])
        )
    elif status == cp_model.UNKNOWN:
        result.status = "TIMEOUT"
    else:
        result.status = "INFEASIBLE"
    print(f"    flips={result.flips}  status={result.status}")
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

ALL_CONFIGS = ("+H1+H2+H3", "+H2+H3", "+H2", "+H1+H2",
               "+H1+H3", "+H3", "+H1", "Naive")


def run_ablation(predictions_file: str,
                 complexes_file: str,
                 go_hierarchy_file: str,
                 taxon_id: str,
                 threshold: float,
                 top_k: int,
                 timeout: int,
                 dry_run: bool = False,
                 configs: Optional[Set[str]] = None) -> List[RunResult]:

    mode_tag = "DRY RUN" if dry_run else "SOLVE"
    print(f"\n{'='*60}")
    print(f"  Heuristic ablation ({mode_tag})  |  taxon {taxon_id}  |"
          f"  threshold={threshold}  top_k={top_k}  timeout={timeout}s")
    print(f"{'='*60}\n")

    # --- Load ---
    print("Loading data...")
    predictions  = load_predictions(predictions_file)
    hierarchy    = load_go_hierarchy(go_hierarchy_file)
    homodimers   = load_homodimer_terms(complexes_file)

    annotated, _ = get_annotated_predictions(predictions, threshold)
    flip_costs   = compute_flip_costs(predictions, threshold)

    all_complexes        = get_heteromeric_complexes(homodimers, hierarchy)
    incoherent_complexes = get_incoherent_complexes(annotated, all_complexes)

    print(f"  Proteins            : {len(predictions):,}")
    print(f"  All complex terms   : {len(all_complexes):,}")
    print(f"  Incoherent complexes: {len(incoherent_complexes):,}")
    print(f"  GO hierarchy edges  : {sum(len(v) for v in hierarchy.values()):,}\n")

    results: List[RunResult] = []

    # Pre-compute participating protein sets for the various configs.
    # top-k for incoherent only
    part_incoh_topk = get_participating_proteins(
        annotated, incoherent_complexes, predictions, top_k=top_k)
    # all proteins for incoherent only
    part_incoh_all = get_participating_proteins(
        annotated, incoherent_complexes, predictions, top_k=None)
    # top-k for all complexes
    part_all_topk = get_participating_proteins(
        annotated, all_complexes, predictions, top_k=top_k)
    # all proteins for all complexes (sparse, no top-k)
    part_all_all = get_participating_proteins(
        annotated, all_complexes, predictions, top_k=None)

    # Pre-compute sparse pair sets.
    pairs_incoh_topk = get_go_protein_pairs(annotated, hierarchy, part_incoh_topk)
    pairs_incoh_all  = get_go_protein_pairs(annotated, hierarchy, part_incoh_all)
    pairs_all_topk   = get_go_protein_pairs(annotated, hierarchy, part_all_topk)
    pairs_all_all    = get_go_protein_pairs(annotated, hierarchy, part_all_all)

    # Run configs from fastest to slowest.
    def _run(label):
        return configs is None or label in configs

    # --- +H1+H2+H3 : incoherent + sparse + top-k ---
    label = "+H1+H2+H3"
    if _run(label):
        print(f"--- {label} (incoherent complexes + sparse vars + top-{top_k}) ---")
        if dry_run:
            r = count_sparse_config(annotated, hierarchy, pairs_incoh_topk,
                                    part_incoh_topk, label)
        else:
            r = solve_sparse_vars(annotated, hierarchy, flip_costs,
                                  part_incoh_topk, pairs_incoh_topk, timeout, label)
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H2+H3 : all complexes + sparse + top-k ---
    label = "+H2+H3"
    if _run(label):
        print(f"--- {label} (all complexes + sparse vars + top-{top_k}) ---")
        if dry_run:
            r = count_sparse_config(annotated, hierarchy, pairs_all_topk,
                                    part_all_topk, label)
        else:
            try:
                r = solve_sparse_vars(annotated, hierarchy, flip_costs,
                                      part_all_topk, pairs_all_topk, timeout, label)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H2 : all complexes + sparse + all proteins ---
    label = "+H2"
    if _run(label):
        print(f"--- {label} (all complexes + sparse vars, all proteins) ---")
        if dry_run:
            r = count_sparse_config(annotated, hierarchy, pairs_all_all,
                                    part_all_all, label)
        else:
            try:
                r = solve_sparse_vars(annotated, hierarchy, flip_costs,
                                      part_all_all, pairs_all_all, timeout, label)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H1+H2 : incoherent + sparse + all proteins ---
    label = "+H1+H2"
    if _run(label):
        print(f"--- {label} (incoherent complexes + sparse vars, all proteins) ---")
        if dry_run:
            r = count_sparse_config(annotated, hierarchy, pairs_incoh_all,
                                    part_incoh_all, label)
        else:
            try:
                r = solve_sparse_vars(annotated, hierarchy, flip_costs,
                                      part_incoh_all, pairs_incoh_all, timeout, label)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H1+H3 : incoherent + full vars + top-k coherence ---
    label = "+H1+H3"
    if _run(label):
        print(f"--- {label} (incoherent complexes + full vars + top-{top_k} coherence) ---")
        if dry_run:
            r = count_full_vars_config(annotated, hierarchy, incoherent_complexes,
                                       label, participating_proteins=part_incoh_topk)
        else:
            try:
                r = solve_full_vars(annotated, hierarchy, flip_costs,
                                    incoherent_complexes, timeout, label,
                                    participating_proteins=part_incoh_topk)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H3 : all complexes + full vars + top-k coherence ---
    label = "+H3"
    if _run(label):
        print(f"--- {label} (all complexes + full vars + top-{top_k} coherence) ---")
        if dry_run:
            r = count_full_vars_config(annotated, hierarchy, all_complexes,
                                       label, participating_proteins=part_all_topk)
        else:
            try:
                r = solve_full_vars(annotated, hierarchy, flip_costs,
                                    all_complexes, timeout, label,
                                    participating_proteins=part_all_topk)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- +H1 : incoherent + full vars ---
    label = "+H1"
    if _run(label):
        print(f"--- {label} (incoherent complexes only, full variable set) ---")
        if dry_run:
            r = count_full_vars_config(annotated, hierarchy, incoherent_complexes,
                                       label)
        else:
            try:
                r = solve_full_vars(annotated, hierarchy, flip_costs,
                                    incoherent_complexes, timeout, label)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    # --- Naive : all complexes + full vars ---
    label = "Naive"
    if _run(label):
        print(f"--- {label} (all complexes, full variable set) ---")
        if dry_run:
            r = count_full_vars_config(annotated, hierarchy, all_complexes, label)
        else:
            try:
                r = solve_full_vars(annotated, hierarchy, flip_costs,
                                    all_complexes, timeout, label)
            except MemoryError:
                r = RunResult(label)
                r.status = "OOM"
                print("    OUT OF MEMORY\n")
        results.append(r)
        print(f"    vars={r.n_vars:,}  constraints={r.n_constraints:,}  status={r.status}\n")

    return results


def print_table(results: List[RunResult], taxon_id: str) -> None:
    print(f"\n{'='*60}")
    print(f"  Summary  |  taxon {taxon_id}")
    print(f"{'='*60}")
    header = (f"{'Config':<15} {'#vars':>12} {'#constraints':>14}"
              f" {'setup(s)':>9} {'solve(s)':>9} {'flips':>7} {'status'}")
    print(header)
    print("-" * len(header))
    for r in results:
        vars_str  = f"{r.n_vars:,}"        if r.n_vars        else "—"
        cons_str  = f"{r.n_constraints:,}"  if r.n_constraints else "—"
        setup_str = f"{r.setup_s:.2f}"      if r.setup_s       else "—"
        solve_str = f"{r.solve_s:.2f}"      if r.solve_s       else "—"
        if r.status in ("OOM", "NOT_RUN", "DRY_RUN"):
            flips_str = "—"
        else:
            flips_str = str(r.flips)
        print(f"{r.config:<15} {vars_str:>12} {cons_str:>14} {setup_str:>9}"
              f" {solve_str:>9} {flips_str:>7}  {r.status}")
    print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Heuristic ablation for complex coherence CP-SAT solver."
    )
    parser.add_argument("--predictions",  required=True,
                        help="Path to predictions TSV file for a single organism.")
    parser.add_argument("--complexes",    required=True,
                        help="Path to protein_complexes.tsv (homodimer classification).")
    parser.add_argument("--go_hierarchy", required=True,
                        help="Path to go_hierarchy_cc.tsv.")
    parser.add_argument("--taxon_id",     default="unknown",
                        help="NCBI taxon ID (used for display only).")
    parser.add_argument("--threshold",    type=float, default=0.3,
                        help="Prediction score threshold (default: 0.3).")
    parser.add_argument("--top_k",        type=int,   default=50,
                        help="Top-k proteins per complex for H3 (default: 50).")
    parser.add_argument("--timeout",      type=int,   default=300,
                        help="CP-SAT solver timeout in seconds (default: 300).")
    parser.add_argument("--dry_run",      action="store_true",
                        help="Only count variables and constraints; skip model "
                             "construction and solving.")
    parser.add_argument("--configs", default=None,
                        help="Comma-separated subset of configs to run "
                             "(e.g. '+H2' or '+H2,+H3'). Default: all 8.")
    args = parser.parse_args()

    configs_set = None
    if args.configs:
        configs_set = {c.strip() for c in args.configs.split(",") if c.strip()}
        unknown = configs_set - set(ALL_CONFIGS)
        if unknown:
            parser.error(f"Unknown config(s): {sorted(unknown)}. "
                         f"Valid: {list(ALL_CONFIGS)}")

    results = run_ablation(
        predictions_file=args.predictions,
        complexes_file=args.complexes,
        go_hierarchy_file=args.go_hierarchy,
        taxon_id=args.taxon_id,
        threshold=args.threshold,
        top_k=args.top_k,
        timeout=args.timeout,
        dry_run=args.dry_run,
        configs=configs_set,
    )
    print_table(results, args.taxon_id)
