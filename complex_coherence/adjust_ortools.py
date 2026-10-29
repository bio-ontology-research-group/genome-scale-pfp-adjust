"""
Adjust the predictions to be complex coherent using min Flip SAT optimizer with OR-Tools.

Inputs:
- complex_predictions.tsv: tab-separated file with protein name, GO term, and prediction score.
  Example: protein_name\tGO:term|score\tGO:term|score...
- protein_complexes.tsv: tab-separated file with GO term, classification, and definition.
- go_hierarchy.tsv: tab-separated file with child GO term and parent GO term.

Outputs:
- complex_adjusted_predictions.tsv: tab-separated file with protein name, GO term, and adjusted prediction score.

This version uses Google OR-Tools CP-SAT.
"""

import csv
from collections import defaultdict
from typing import Dict, Set, Tuple, List
from ortools.sat.python import cp_model
import time
import argparse

# Root terms for each GO namespace (should be excluded from the problem)
ROOT_TERMS = {'GO:0005575', 'GO:0008150', 'GO:0003674'}  # CC, BP, MF


def load_predictions(file_path):
    """
    Load the predictions from the file.
    """
    predictions = defaultdict(dict)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if not row or not row[0]:
                continue
            protein_name = row[0]

            predictions[protein_name] = {
                go_term: float(score) for go_term_score in row[1:] for go_term, score in [go_term_score.split('|')]
            }
    return predictions


def load_homodimer_terms(file_path):
    """
    Load the homodimer complexes terms from the file.
    Expected format: GO_ID\tClassification\tDefinition

    returns a set of GO terms that are homodimeric complexes.
    """
    homodimer_terms = set()
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            go_term = row[0]
            if go_term == "GO_term":
                continue  # skip header row
            classification = row[1]

            if classification in ['h', 'a']: # 'h' for homodimer, 'a' for ambiguous I think. Following GAEF's complex_classifier.py.
                homodimer_terms.add(go_term)
    return homodimer_terms

# def load_heteromeric_complexes(file_path):
#     """
#     Load the heteromeric complexes from the file.
#     Expected format: GO_ID\tClassification\tDefinition

#     returns a set of GO terms that are heteromeric complexes.
#     """
#     heteromeric_complexes = set()
#     with open(file_path, 'r') as file:
#         reader = csv.reader(file, delimiter='\t')
#         for row in reader:
#             go_term = row[0]
#             if go_term == "GO_term":
#                 continue  # skip header row
#             classification = row[1]
#             definition = row[2]
#             if classification not in ['h', 'a']:
#                 # would be 'n'. TODO: follow up on this.
#                 heteromeric_complexes.add(go_term)
#     return heteromeric_complexes


# def expand_heteromeric_complexes(heteromeric_complexes: Set[str], go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
#     """
#     Expand the heteromeric complexes to include all children of the complexes.
#     """
#     expanded_heteromeric_complexes = set()
#     for heteromeric_complex in heteromeric_complexes:
#         all_children = get_all_children(heteromeric_complex, go_hierarchy)
#         expanded_heteromeric_complexes.update(all_children)

#     return expanded_heteromeric_complexes.union(heteromeric_complexes)

def get_heteromeric_complexes(homodimer_terms: Set[str], annotated_go_terms: Set[str], go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    """
    Get the heteromeric complexes from the homodimer terms and annotated go terms.
    """
    MACROMOLECULAR_COMPLEX = "GO:0032991"
    all_complex_terms = get_all_children(MACROMOLECULAR_COMPLEX, go_hierarchy)
    all_complex_terms.add(MACROMOLECULAR_COMPLEX)

    print(f"[INFO] all_complex_terms: {len(all_complex_terms)}")

    heteromeric_complexes = all_complex_terms.difference(homodimer_terms)

    return heteromeric_complexes


def get_all_children(go_term: str, go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    return get_descendants(go_term, go_hierarchy)


def load_go_hierarchy(file_path):
    """
    Load the GO hierarchy from the file.
    """
    hierarchy = defaultdict(set)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            child_go_term = row[0]
            parent_go_term = row[1]
            hierarchy[child_go_term].add(parent_go_term)
    return hierarchy


def compute_flip_cost(predictions: Dict[str, Dict[str, float]],
                       fold_threshold: float = 0.5) -> Tuple[Dict[str, Dict[str, float]], int]:
    """
    Compute flip cost for the predictions.
    """

    # Create nested defaultdict with fold_threshold as default value
    flip_cost = defaultdict(lambda: defaultdict(lambda: fold_threshold))
    num_annotated_predictions = 0
    for protein in predictions:
        for go_term, score in predictions[protein].items():
            flip_cost[protein][go_term] = abs(score - fold_threshold)
        num_annotated_predictions += len(predictions[protein])

    print(f"Computed {num_annotated_predictions} flip costs for {len(flip_cost)} proteins")
    return flip_cost, num_annotated_predictions


def get_annotated_predictions(predictions: Dict[str, Dict[str, float]], fold_threshold: float = 0.5) -> Tuple[Dict[str, Set[str]], Set[str]]:
    """
    Get the annotated predictions from the predictions.
    """
    annotated_predictions = defaultdict(set, {
        protein: {go_term for go_term, score in go_term_scores.items() if score > fold_threshold}
        for protein, go_term_scores in predictions.items()
    })

    annotated_go_terms = set().union(*annotated_predictions.values())
    return annotated_predictions, annotated_go_terms



def solve_sat_ortools_hierarchy(annotated_predictions: Dict[str, Set[str]],
                go_hierarchy: Dict[str, Set[str]] = None,
                flip_costs: Dict[str, Dict[str, float]] = None,
                heteromeric_complexes: Set[str] = None,
                demotion_only: bool = False,
               ) -> Dict[str, Dict[str, bool]]:
    """
    Solve the SAT problem using OR-Tools CP-SAT.

    Args:
        annotated_predictions: Dict of sets of strings, where annotated_predictions[protein] is the set of go terms annotated to protein.
        go_hierarchy: Dict of sets of strings, where go_hierarchy[go_term] is the set of parent go terms of go_term.
        flip_costs: Dict of floats, where flip_costs[protein][go_term] is the cost of flipping the prediction for go_term.
        heteromeric_complexes: Set of strings, the GO terms that are heteromeric complexes.

    Returns:
        Dict of dicts of booleans, where adjusted_predictions[protein][go_term] is True if go_term is flipped for protein.
    """
    setup_start_time = time.time()
    model = cp_model.CpModel()
    
    # Handle None go_hierarchy
    if go_hierarchy is None:
        go_hierarchy = defaultdict(set)
    
    # Handle None flip_costs
    if flip_costs is None:
        flip_costs = defaultdict(lambda: defaultdict(lambda: 1.0))

    heteromeric_complexes = heteromeric_complexes or set()
    # Add GO terms variables
    hierarchy_go_terms = set()
    all_go_terms = set()
    go_term_variables = {}
    # get all go terms from go_hierarchy and annotated_predictions
    for protein in sorted(annotated_predictions.keys()):
        for child_go_term in sorted(go_hierarchy.keys()):
            for parent_go_term in sorted(go_hierarchy.get(child_go_term, ())):
                hierarchy_go_terms.add(child_go_term)
                hierarchy_go_terms.add(parent_go_term)
                all_go_terms.add((child_go_term, protein))
                all_go_terms.add((parent_go_term, protein))

    hierarchy_go_terms.update(set().union(*annotated_predictions.values()))
    hierarchy_go_terms.update(heteromeric_complexes)
    all_go_terms.update((g, p) for g in hierarchy_go_terms for p in annotated_predictions)

    for go_term, protein in sorted(all_go_terms):
        go_term_variables[(go_term, protein)] = model.NewBoolVar(f'go_{go_term}_{protein}')

    print(f"[INFO] Added {len(go_term_variables)} GO term variables")

    # Add GO hierarchy constraints
    for child_go_term, protein in sorted(all_go_terms):
        for parent_go_term in sorted(go_hierarchy.get(child_go_term, ())):
            if (parent_go_term, protein) in all_go_terms:
                # child implies parent: child <= parent
                model.AddImplication(go_term_variables[(child_go_term, protein)], 
                                    go_term_variables[(parent_go_term, protein)])


    # Add heteromeric complex constraints
    for heteromeric_complex in sorted(heteromeric_complexes):
        for protein in sorted(annotated_predictions.keys()):
            if (heteromeric_complex, protein) in go_term_variables:
                other_proteins = sorted([p for p in annotated_predictions.keys() if p != protein])
                other_vars = [go_term_variables[(heteromeric_complex, other_protein)] 
                             for other_protein in other_proteins 
                             if (heteromeric_complex, other_protein) in go_term_variables]
                if other_vars:
                    # If this protein has the complex, at least one other must have it
                    model.AddBoolOr(other_vars).OnlyEnforceIf(go_term_variables[(heteromeric_complex, protein)])
                else:
                    model.Add(go_term_variables[(heteromeric_complex, protein)] == 0)

    print(f"[INFO] Added {len(heteromeric_complexes)} heteromeric complex constraints")


    # Add flip variables
    flip_variables = {}
    flip_go_terms = set()
    for go_term, protein in sorted(all_go_terms):
        if (go_term, protein) not in flip_go_terms:
            flip_go_terms.add((go_term, protein))
            flip_variables[(go_term, protein)] = model.NewBoolVar(f'flip_{go_term}_{protein}')

    print(f"[INFO] Added {len(flip_variables)} flip variables")

    # Add theta and eta constraints
    # Z3 formulation: theta = go_var XOR flip_var (for annotated terms)
    #                 eta = NOT(go_var) XOR flip_var (for non-annotated terms)
    # All theta and eta must be True (enforced via protein -> AND(thetas, etas) and genome -> AND(proteins))
    
    # Since genome and all proteins are forced True, we directly enforce:
    # For annotated terms: go_var XOR flip_var = True  =>  go_var + flip_var == 1
    # For non-annotated terms: NOT(go_var) XOR flip_var = True  =>  go_var == flip_var
    
    for protein in sorted(annotated_predictions.keys()):
        # For annotated GO terms: XOR must be True
        for go_term in sorted(annotated_predictions[protein]):
            if (go_term, protein) in all_go_terms:
                # go_var XOR flip_var = True means exactly one is true
                model.Add(go_term_variables[(go_term, protein)] + flip_variables[(go_term, protein)] == 1)
        
        # For non-annotated GO terms: NOT(go_var) XOR flip_var = True
        for go_term in sorted(hierarchy_go_terms):
            if go_term not in annotated_predictions[protein] and (go_term, protein) in all_go_terms:
                # NOT(go_var) XOR flip_var = True
                # This means: (1 - go_var) XOR flip_var = True
                # Which is: (1 - go_var) + flip_var = 1 (one of them is true)
                # Simplifies to: go_var = flip_var
                model.Add(go_term_variables[(go_term, protein)] == flip_variables[(go_term, protein)])
                if demotion_only:
                    model.Add(go_term_variables[(go_term, protein)] == 0)

    print(f"[INFO] Added {len(annotated_predictions)} protein variables")

    # Add cost constraints - scale to integers with high precision
    # Use 1000000 scale factor to preserve precision down to 0.000001
    cost_scale = 1000000
    cost_terms = []
    for go_term, protein in sorted(all_go_terms):
        if (go_term, protein) in flip_variables:
            cost = max(1, int(flip_costs[protein][go_term] * cost_scale))
            # Ensure minimum cost of 1 to prevent "free" flips
            if cost == 0 and flip_costs[protein][go_term] > 0:
                cost = 1
            # Add cost * flip_var to objective
            cost_terms.append(flip_variables[(go_term, protein)] * cost)
    
    model.Minimize(sum(cost_terms))

    solve_start_time = time.time()
    print(f"[INFO] Setup time: {solve_start_time - setup_start_time:.2f}s")
    print("[INFO] Solving SAT problem...")
    
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 42
    status = solver.Solve(model)
    solve_time = time.time() - solve_start_time


    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adjusted_predictions = defaultdict(dict)
        total_flips = 0
        print(f"[INFO] SAT problem solved successfully in {solve_time:.2f}s")
        for go_term, protein in all_go_terms:
            if (go_term, protein) in flip_variables:
                flip_result = solver.BooleanValue(flip_variables[(go_term, protein)])
                go_term_result = solver.BooleanValue(go_term_variables[(go_term, protein)])

                if flip_result:
                    total_flips += 1
                
                adjusted_predictions[protein][go_term] = go_term_result
        print(f"[INFO]   Total flips: {total_flips}")

        opt_cost = solver.ObjectiveValue() / cost_scale
        print(f"[INFO]   Optimal cost: {opt_cost}")

        return adjusted_predictions, total_flips
    else:
        raise RuntimeError(f"CP-SAT did not produce a valid repair: {solver.StatusName(status)}")


def get_ancestors(go_term: str, go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    ancestors = set()
    pending = list(go_hierarchy.get(go_term, ()))
    while pending:
        term = pending.pop()
        if term not in ancestors:
            ancestors.add(term)
            pending.extend(go_hierarchy.get(term, ()))
    ancestors.discard(go_term)
    return ancestors


def get_descendants(go_term: str, go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    # Traverse each reachable node once; recursive DAG expansion revisits shared
    # descendants exponentially on the real GO graph.
    children = defaultdict(set)
    for child, parents in go_hierarchy.items():
        for parent in parents:
            children[parent].add(child)
    descendants = set()
    pending = list(children.get(go_term, ()))
    while pending:
        term = pending.pop()
        if term not in descendants:
            descendants.add(term)
            pending.extend(children.get(term, ()))
    descendants.discard(go_term)
    return descendants


def get_go_protein_pairs(
    annotated_predictions: Dict[str, Set[str]], 
    go_hierarchy: Dict[str, Set[str]],
    participating_proteins: Dict[str, List[str]],
) -> Set[Tuple[str, str]]:
    """
    Get the (go_term, protein) pairs to create their variables in the SAT problem.
    Conceptually, from participating proteins, a predicted complex term is consdiered with all its descendants (for potential demotion), and a non-predicted complex term is considered with all its ancestors (for potential promotion) from the set of originally predicted go terms.
    """
    go_protein_pairs = set()

    for complex_go_term in participating_proteins:
        desc_terms = get_descendants(complex_go_term, go_hierarchy)
        anc_terms = get_ancestors(complex_go_term, go_hierarchy)
        for protein in participating_proteins[complex_go_term]:
            go_protein_pairs.add((complex_go_term, protein))
            if complex_go_term in annotated_predictions[protein]:
                for descendant_go_term in desc_terms:
                    if descendant_go_term in annotated_predictions[protein]:
                        go_protein_pairs.add((descendant_go_term, protein))
            else:
                for parent_go_term in anc_terms:
                    if parent_go_term not in annotated_predictions[protein]:
                        go_protein_pairs.add((parent_go_term, protein))
    return go_protein_pairs


def prepare_sparse_domain(annotated_predictions, go_hierarchy, participating_proteins, complex_terms):
    """Close the repair domain over hierarchy and affected complex constraints.

    Terms outside the model remain fixed. If changing a previously coherent
    complex can leave one fixed member, include that member as well. This avoids
    creating a new singleton while repairing an existing one.
    """
    pairs = get_go_protein_pairs(annotated_predictions, go_hierarchy, participating_proteins)
    ancestors, descendants = {}, {}
    members = defaultdict(set)
    for protein, terms in annotated_predictions.items():
        for term in terms & complex_terms:
            members[term].add(protein)
        for term in terms:
            if any(parent not in terms for parent in go_hierarchy.get(term, ())):
                pairs.add((term, protein))
    pending = list(pairs)
    while pending:
        term, protein = pending.pop()
        ancestors.setdefault(term, None)
        if ancestors[term] is None:
            ancestors[term] = get_ancestors(term, go_hierarchy)
        additions = {(parent, protein) for parent in ancestors[term]
                     if parent not in annotated_predictions[protein]}
        if term in annotated_predictions[protein]:
            if term not in descendants:
                descendants[term] = get_descendants(term, go_hierarchy)
            additions.update((child, protein) for child in descendants[term]
                             if child in annotated_predictions[protein])
        if term in complex_terms:
            fixed = {p for p in members[term] if (term, p) not in pairs}
            if len(fixed) < 2:
                additions.update((term, p) for p in fixed)
        new = additions - pairs
        pairs.update(new)
        pending.extend(sorted(new))
    affected = defaultdict(list)
    for term, protein in sorted(pairs):
        if term in complex_terms:
            affected[term].append(protein)
    return pairs, affected


def solve_sat_ortools_optimized(annotated_predictions: Dict[str, Set[str]],
                go_protein_pairs: Set[Tuple[str, str]],
                go_hierarchy: Dict[str, Set[str]],
                flip_costs: Dict[str, Dict[str, float]],
                participating_proteins: Dict[str, List[str]] = None,
                demotion_only: bool = False,
                violation_penalty: float = None,
                time_limit: float = None,
                solver_stats: dict = None,
                allowed_promotions: set = None,
                allowed_demotions: set = None,
                protected_pairs: set = None,
               ) -> Dict[str, Dict[str, bool]]:
    """
    Solve the SAT problem using OR-Tools CP-SAT with optimized constraints.

    Args:
        annotated_predictions: Dict of sets of strings, where annotated_predictions[protein] is the set of go terms annotated to protein.
        go_protein_pairs: Set of tuples of strings, the (go_term, protein) pairs to create their variables in the SAT problem.
        go_hierarchy: Dict of sets of strings, where go_hierarchy[go_term] is the set of parent go terms of go_term.
        flip_costs: Dict of floats, where flip_costs[protein][go_term] is the cost of flipping the prediction for go_term.
        participating_proteins: Dict of lists of strings, where participating_proteins[complex_go_term] is the list of proteins that participate in the heteromeric complex complex_go_term. (can be top-k or above threshold)
    Returns:
        Dict of dicts of booleans, where adjusted_predictions[protein][go_term] is True if go_term is flipped for protein.
    """
    for term, protein in protected_pairs or ():
        if term not in annotated_predictions.get(protein, set()):
            raise ValueError('Protected pairs must be active in the input')
    setup_start_time = time.time()
    model = cp_model.CpModel()
    
    annotated_go_terms = set()
    for protein in annotated_predictions:
        for go_term in annotated_predictions[protein]:
            annotated_go_terms.add(go_term)

    # Add GO term variables and flip variables
    all_go_terms = set()
    go_term_variables = {}
    flip_variables = {}
    proteins_added = set()
    for go_term, protein in sorted(go_protein_pairs):
        all_go_terms.add((go_term, protein))
        go_term_variables[(go_term, protein)] = model.NewBoolVar(f'go_{go_term}_{protein}')
        # A flip is the complement of an active input, or the decision itself
        # for an inactive input. Using that literal directly avoids a second
        # variable and equality per pair without changing the objective.
        variable = go_term_variables[(go_term, protein)]
        flip_variables[(go_term, protein)] = variable.Not() if go_term in annotated_predictions[protein] else variable
        proteins_added.add(protein)

    print(f"[INFO] Added {len(go_term_variables)} decision variables with derived flip literals")
    print(f"[INFO] Added {len(proteins_added)} proteins")

    all_go_terms = sorted(all_go_terms) # sort for deterministic output

    # Add GO hierarchy constraints
    for child_go_term, protein in all_go_terms:
        for parent_go_term in sorted(go_hierarchy.get(child_go_term, ())):
            if (parent_go_term, protein) in go_term_variables:
                # child implies parent: child <= parent
                model.AddImplication(go_term_variables[(child_go_term, protein)], 
                                    go_term_variables[(parent_go_term, protein)])
            elif parent_go_term not in annotated_predictions[protein]:
                model.Add(go_term_variables[(child_go_term, protein)] == 0)
    # Boundary terms outside the sparse model retain their original state.
    children = defaultdict(set)
    for child, parents in go_hierarchy.items():
        for parent in parents:
            children[parent].add(child)
    for parent, protein in all_go_terms:
        if any(child in annotated_predictions[protein] and (child, protein) not in go_term_variables
               for child in children[parent]):
            model.Add(go_term_variables[(parent, protein)] == 1)


    # Add heteromeric complex constraints
    violation_variables = []
    for heteromeric_complex in sorted(participating_proteins.keys()):
        complex_vars = [go_term_variables[(heteromeric_complex, protein)] 
                        for protein in participating_proteins[heteromeric_complex]
                        if (heteromeric_complex, protein) in go_term_variables]
        fixed_members = sum(heteromeric_complex in terms and (heteromeric_complex, protein) not in go_term_variables
                            for protein, terms in annotated_predictions.items())
        if complex_vars:
            if violation_penalty is not None:
                singleton = model.NewBoolVar(f'singleton_{heteromeric_complex}')
                model.Add(sum(complex_vars) + fixed_members == 1).OnlyEnforceIf(singleton)
                model.Add(sum(complex_vars) + fixed_members != 1).OnlyEnforceIf(singleton.Not())
                violation_variables.append(singleton)
                continue
            no_members = model.NewBoolVar(f'complex_{heteromeric_complex}_has_zero_members')
            model.Add(sum(complex_vars) + fixed_members == 0).OnlyEnforceIf(no_members)
            model.Add(sum(complex_vars) + fixed_members >= 2).OnlyEnforceIf(no_members.Not())

    print(f"[INFO] Added {len(participating_proteins)} heteromeric complex constraints")

    # Flip literals already encode whether each input decision changed.
    for go_term, protein in all_go_terms:
        if demotion_only and go_term not in annotated_predictions[protein]:
            model.Add(go_term_variables[(go_term, protein)] == 0)
        pair = (go_term, protein)
        active = go_term in annotated_predictions[protein]
        if not active and allowed_promotions is not None and pair not in allowed_promotions:
            model.Add(go_term_variables[pair] == 0)
        if active and allowed_demotions is not None and pair not in allowed_demotions:
            model.Add(go_term_variables[pair] == 1)
        if pair in (protected_pairs or ()):
            if not active:
                raise ValueError('Protected pairs must be active in the input')
            model.Add(go_term_variables[pair] == 1)

    # Add cost constraints - scale to integers with high precision
    # Use 1000000 scale factor to preserve precision down to 0.000001
    cost_scale = 1000000
    cost_terms = []
    for go_term, protein in all_go_terms:
        if (go_term, protein) in flip_variables:
            cost = max(1, int(flip_costs[protein][go_term] * cost_scale))
            cost_terms.append(flip_variables[(go_term, protein)] * cost)
    if violation_penalty is not None:
        if violation_penalty < 0:
            raise ValueError('violation_penalty must be nonnegative')
        cost_terms.extend(var * round(violation_penalty * cost_scale) for var in violation_variables)

    # A feasible demotion repair is a useful starting point for symmetric search.
    # It does not change the objective or constrain the optimized solution.
    hinted = {(g, p): g in annotated_predictions[p] for g, p in all_go_terms}
    changed = True
    while changed:
        changed = False
        for (g, p), present in list(hinted.items()):
            if present and any(not hinted.get((parent, p), parent in annotated_predictions[p])
                               for parent in go_hierarchy.get(g, ())):
                hinted[g, p] = False
                changed = True
        for term, proteins in participating_proteins.items():
            mutable_members = [(term, p) for p in proteins if hinted.get((term, p), False)]
            fixed_count = sum(term in terms and (term, p) not in hinted
                              for p, terms in annotated_predictions.items())
            if fixed_count + len(mutable_members) == 1 and mutable_members:
                hinted[mutable_members[0]] = False
                changed = True
    if allowed_promotions is not None or allowed_demotions is not None:
        # The ancestor-closed input is feasible for the abstaining formulation;
        # the ordinary demotion hint may conflict with protected annotations.
        hinted = {(g, p): g in annotated_predictions[p] for g, p in all_go_terms}
    for pair, present in hinted.items():
        model.AddHint(go_term_variables[pair], int(present))
    
    model.Minimize(sum(cost_terms))

    solve_start_time = time.time()
    print(f"[INFO] Setup time: {solve_start_time - setup_start_time:.2f}s")
    print("[INFO] Solving SAT problem...")
    
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    if time_limit is not None:
        solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.random_seed = 42
    status = solver.Solve(model)
    solve_time = time.time() - solve_start_time
    if solver_stats is not None:
        solver_stats.update(status=solver.StatusName(status), seconds=solve_time,
                            objective=solver.ObjectiveValue()/cost_scale,
                            bound=solver.BestObjectiveBound()/cost_scale,
                            optimal=(status == cp_model.OPTIMAL))


    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adjusted_predictions = defaultdict(dict)
        total_flips = 0
        print(f"[INFO] SAT problem solved successfully in {solve_time:.2f}s")
        for go_term, protein in all_go_terms:
            if (go_term, protein) in flip_variables:
                flip_result = solver.BooleanValue(flip_variables[(go_term, protein)])
                go_term_result = solver.BooleanValue(go_term_variables[(go_term, protein)])

                if flip_result:
                    total_flips += 1
                
                adjusted_predictions[protein][go_term] = go_term_result
        print(f"[INFO]   Total flips: {total_flips}")

        opt_cost = solver.ObjectiveValue() / cost_scale
        print(f"[INFO]   Best cost: {opt_cost}; status={solver.StatusName(status)}")

        return adjusted_predictions, total_flips
    else:
        raise RuntimeError(f"CP-SAT did not produce a valid repair: {solver.StatusName(status)}")


def get_participating_proteins(annotated_predictions: Dict[str, Set[str]], complex_terms: Set[str], predictions: Dict[str, Dict[str, float]], top_k = None, participating_threshold = None, demotion_only: bool = False) -> Dict[str, List[str]]:
    """
    Get the participating proteins for each heteromeric complex.
    If top_k value is provided, consider only incoherent complexes and return the top-k proteins for each complex.
    If participating_threshold value is provided, return the proteins above the threshold for each complex.
    If neither are provided, return all proteins for each complex.
    """
    participating_proteins = defaultdict(list)
    if demotion_only:
        for complex_go_term in complex_terms:
            proteins_with_complex = [p for p in annotated_predictions.keys()
                                     if complex_go_term in annotated_predictions[p]]
            if len(proteins_with_complex) != 1:
                continue
            participating_proteins[complex_go_term] = proteins_with_complex
    elif top_k is not None:
        for complex_go_term in complex_terms:
            # count how many proteins have the complex term
            proteins_with_complex = [p for p in annotated_predictions.keys() 
                                     if complex_go_term in annotated_predictions[p]]
            if len(proteins_with_complex) != 1:
                continue # only consider incoherent complexes
            
            complex_predictions_scores = [(p, predictions[p].get(complex_go_term, 0.0)) for p in predictions.keys()]
            sorted_proteins = sorted(complex_predictions_scores, key=lambda x: (-x[1], x[0]))
            participating_proteins[complex_go_term] = [p for p, score in sorted_proteins[:top_k]]

            
    elif participating_threshold is not None:
        for complex_go_term in sorted(complex_terms):
            active_members = {p for p, terms in annotated_predictions.items() if complex_go_term in terms}
            if len(active_members) != 1:
                continue
            candidates = {p for p, scores in predictions.items() if scores.get(complex_go_term, 0) > participating_threshold}
            participating_proteins[complex_go_term] = sorted(active_members | candidates)
    else:
        for complex_go_term in complex_terms:
            # count how many proteins have the complex term
            proteins_with_complex = [p for p in annotated_predictions.keys() 
                                     if complex_go_term in annotated_predictions[p]]
            if len(proteins_with_complex) != 1:
                continue # only consider incoherent complexes

            participating_proteins[complex_go_term] = [protein for protein in annotated_predictions.keys()]
    return participating_proteins


def get_final_predictions(adjusted_predictions: Dict[str, Dict[str, bool]], predictions: Dict[str, Dict[str, float]], annotated_go_terms: Set[str], fold_threshold: float = 0.5, epsilon: float = 0.001) -> Dict[str, Dict[str, float]]:
    """
    Export Boolean repairs while preserving every unmodified score.

    Presence means score > fold_threshold. Demotions are capped at the
    threshold; promotions are raised by epsilon. annotated_go_terms is retained
    for API compatibility, but must not filter the original score vocabulary.
    """
    if not 0 <= fold_threshold < 1 or epsilon <= 0:
        raise ValueError('Require threshold in [0,1) and positive epsilon')
    final_predictions = defaultdict(dict)
    for protein, p_scores in predictions.items():
        adj_p = adjusted_predictions.get(protein)
        if not adj_p:
            final_predictions[protein] = {
                go_term: score
                for go_term, score in p_scores.items()
            }
        else:
            out_p = {
                go_term: score
                for go_term, score in p_scores.items()
            }
            for go_term, is_present in adj_p.items():
                orig_score = p_scores.get(go_term, 0.0)
                if bool(is_present) == (orig_score > fold_threshold):
                    continue
                if is_present:
                    final_score = max(orig_score, min(1.0, fold_threshold + epsilon))
                else:
                    final_score = min(orig_score, fold_threshold)
                if final_score > 0:
                    out_p[go_term] = final_score
                else:
                    out_p.pop(go_term, None)
            final_predictions[protein] = out_p

    return final_predictions


def save_predictions(predictions: Dict[str, Dict[str, float]], file_path: str):
    """
    Save the adjusted predictions to a file.
    """
    with open(file_path, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        for protein, go_terms in predictions.items():
            go_terms_scores = [f"{go_term}|{score:.17g}" for go_term, score in go_terms.items()]

            writer.writerow([protein] + go_terms_scores)



def main(
    predictions_file: str,
    complexes_file: str,
    go_hierarchy_file: str,
    output_file: str,
    threshold: float = 0.5,
    optimized: bool = False,
    top_k: int = None,
    participating_threshold: float = None,
    demotion_only: bool = False,
):
    """
    Main function to run the complex coherence adjustment script.
    """
    
    if not 0 <= threshold < 1:
        raise ValueError('threshold must be in [0, 1)')
    if top_k is not None and top_k < 1:
        raise ValueError('top_k must be positive')
    print("Loading data...")
    predictions = load_predictions(predictions_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file)
    homodimer_terms = load_homodimer_terms(complexes_file)
    # complexes = load_heteromeric_complexes(complexes_file)
    # expanded_complexes = expand_heteromeric_complexes(complexes, go_hierarchy)
    print(f"Loaded {len(predictions)} proteins")
    print(f"Loaded {len(homodimer_terms)} homodimer terms")
    print(f"Loaded {len(go_hierarchy)} GO hierarchy relationships")

    annotated_predictions, annotated_go_terms = get_annotated_predictions(predictions, threshold)

    heteromeric_complexes = get_heteromeric_complexes(homodimer_terms, annotated_go_terms, go_hierarchy)
    print(f"[INFO] Number of Heteromeric Complexes: {len(heteromeric_complexes)}")


    print("Setting up SAT problem...")
    if optimized:
        participating_proteins = get_participating_proteins(annotated_predictions, heteromeric_complexes, predictions, top_k=top_k, participating_threshold=participating_threshold, demotion_only=demotion_only)

        print(f"[INFO] Number of Incoherent Complexes: {len(participating_proteins)}")
        go_protein_pairs, participating_proteins = prepare_sparse_domain(
            annotated_predictions, go_hierarchy, participating_proteins, heteromeric_complexes)

        if not go_protein_pairs:
            print("[INFO] No participating proteins found, returning original predictions")
            if output_file != predictions_file:
                save_predictions(predictions, output_file)
            else:
                print("[INFO] Output file is the same as the predictions file, not saving")
            print("Done!")
            return 0, sum(map(len, annotated_predictions.values()))

        # compute after checking if there are participating proteins
        flip_costs, num_annotated_predictions = compute_flip_cost(predictions, threshold)
        
        adjusted_predictions, total_flips = solve_sat_ortools_optimized(annotated_predictions, 
        go_protein_pairs,
        go_hierarchy, 
        flip_costs,
        participating_proteins,
        demotion_only=demotion_only)
    else:
        flip_costs, num_annotated_predictions = compute_flip_cost(predictions, threshold)

        adjusted_predictions, total_flips = solve_sat_ortools_hierarchy(annotated_predictions, 
        go_hierarchy, 
        flip_costs,
        heteromeric_complexes, demotion_only=demotion_only)

    num_annotated_predictions = sum(map(len, annotated_predictions.values()))
    print(f"[INFO] flips per active input prediction: {total_flips / num_annotated_predictions * 100 if num_annotated_predictions else 0:.2f}%")

    final_predictions = get_final_predictions(adjusted_predictions, predictions, annotated_go_terms, threshold)
    
    print(f"\nSaving final predictions to {output_file}...")
    save_predictions(final_predictions, output_file)
    print("Done!")

    return total_flips, num_annotated_predictions


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Adjust the predictions to be complex coherent using min Flip SAT optimizer with OR-Tools.')
    parser.add_argument('--predictions', required=True, help='Path to complex_predictions.tsv')
    parser.add_argument('--complexes', required=True, help='Path to protein_complexes.tsv')
    parser.add_argument('--go_hierarchy', default="data/go_hierarchy_cc.tsv", help='Path to go_hierarchy.tsv')
    parser.add_argument('--output', required=True, help='Path to output complex_adjusted_predictions.tsv')
    parser.add_argument('--threshold', type=float, default=0.5, help='Prediction threshold (default: 0.5)')
    parser.add_argument('--optimized', action='store_true', help='Use optimized constraints')
    parser.add_argument('--top_k', type=int, default=None, help='Top-k proteins to participate in each complex (default: None)')
    parser.add_argument('--participating_threshold', type=float, default=None, help='Participating threshold for proteins to participate in each complex (default: None)')
    parser.add_argument('--demotion-only', '--demotion_only', dest='demotion_only', action='store_true', help='Forbid 0->1 promotions; resolve incoherent complexes by demotion only')
    args = parser.parse_args()
    main(
        predictions_file=args.predictions,
        complexes_file=args.complexes,
        go_hierarchy_file=args.go_hierarchy,
        output_file=args.output,
        threshold=args.threshold,
        optimized=args.optimized,
        top_k=args.top_k,
        participating_threshold=args.participating_threshold,
        demotion_only=args.demotion_only,
    )
