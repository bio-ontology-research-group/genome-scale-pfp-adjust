"""
Adjust the predictions to be complex coherent using min Flip SAT optimizer with OR-Tools.

Inputs:
- complex_predictions.tsv: tab-separated file with protein name, GO term, and prediction score.
  Example: protein_name\tGO:term|score\tGO:term|score...
- protein_complexes.tsv: tab-separated file with GO term, classification, and definition.
- go_hierarchy.tsv: tab-separated file with child GO term and parent GO term.

Outputs:
- complex_adjusted_predictions.tsv: tab-separated file with protein name, GO term, and adjusted prediction score.

This version uses Google OR-Tools CP-SAT for 10-50x speedup over Z3.
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
            if len(row) < 2:
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
    """
    Get all children of the go term using a breadth-first search.
    """
    children = set()
    queue = [go_term]
    while queue:
        current_go_term = queue.pop(0)
        # Sort to ensure deterministic order
        sorted_children = sorted([child for child in go_hierarchy if current_go_term in go_hierarchy[child]])
        queue.extend([child for child in sorted_children if child not in children])
        children.update(sorted_children)
    return children


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
        flip_costs = defaultdict(dict)

    # Add GO terms variables
    hierarchy_go_terms = set()
    all_go_terms = set()
    go_term_variables = {}
    # get all go terms from go_hierarchy and annotated_predictions
    for protein in sorted(annotated_predictions.keys()):
        for child_go_term in sorted(go_hierarchy.keys()):
            for parent_go_term in sorted(go_hierarchy[child_go_term]):
                hierarchy_go_terms.add(child_go_term)
                hierarchy_go_terms.add(parent_go_term)
                all_go_terms.add((child_go_term, protein))
                all_go_terms.add((parent_go_term, protein))

    for go_term, protein in sorted(all_go_terms):
        go_term_variables[(go_term, protein)] = model.NewBoolVar(f'go_{go_term}_{protein}')

    print(f"[INFO] Added {len(go_term_variables)} GO term variables")

    # Add GO hierarchy constraints
    for child_go_term, protein in sorted(all_go_terms):
        for parent_go_term in sorted(go_hierarchy[child_go_term]):
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

    print(f"[INFO] Added {len(annotated_predictions)} protein variables")

    # Add cost constraints - scale to integers with high precision
    # Use 1000000 scale factor to preserve precision down to 0.000001
    cost_scale = 1000000
    cost_terms = []
    for go_term, protein in sorted(all_go_terms):
        if (go_term, protein) in flip_variables:
            cost = int(flip_costs[protein][go_term] * cost_scale)
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
        print(f"[INFO] SAT problem is unsatisfiable after {solve_time:.2f}s")
        print("[INFO] Returning original predictions")
        # Return original annotations as dict of bools
        result = defaultdict(dict)
        for protein, terms in annotated_predictions.items():
            for term in terms:
                result[protein][term] = True
        return result, 0


def get_ancestors(go_term: str, go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    """
    Get the ancestors of a go term.
    
    Args:
        go_term: The GO term to get ancestors for
        go_hierarchy: Dictionary mapping child GO terms to sets of parent GO terms
    
    Returns:
        Set of ancestor GO terms
    """
    if go_term in ROOT_TERMS:
        return set()
    ancestors = set()
    for parent_go_term in go_hierarchy[go_term]:
        ancestors.add(parent_go_term)
        ancestors.update(get_ancestors(parent_go_term, go_hierarchy))
    return ancestors

def get_descendants(go_term: str, go_hierarchy: Dict[str, Set[str]]) -> Set[str]:
    """
    Get the descendants of a go term.
    """
    descendants = set()
    for child_go_term in [child for child in go_hierarchy if go_term in go_hierarchy[child]]:
        descendants.add(child_go_term)
        descendants.update(get_descendants(child_go_term, go_hierarchy))
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
        for protein in participating_proteins[complex_go_term]:
            go_protein_pairs.add((complex_go_term, protein))
            if complex_go_term in annotated_predictions[protein]:
                for descendant_go_term in get_descendants(complex_go_term, go_hierarchy):
                    if descendant_go_term in annotated_predictions[protein]:
                        go_protein_pairs.add((descendant_go_term, protein))
            else:
                for parent_go_term in get_ancestors(complex_go_term, go_hierarchy):
                    if parent_go_term not in annotated_predictions[protein]:
                        go_protein_pairs.add((parent_go_term, protein))
    return go_protein_pairs


def solve_sat_ortools_optimized(annotated_predictions: Dict[str, Set[str]],
                go_protein_pairs: Set[Tuple[str, str]],
                go_hierarchy: Dict[str, Set[str]],
                flip_costs: Dict[str, Dict[str, float]],
                participating_proteins: Dict[str, List[str]] = None,
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
    for go_term, protein in go_protein_pairs:
        all_go_terms.add((go_term, protein))
        go_term_variables[(go_term, protein)] = model.NewBoolVar(f'go_{go_term}_{protein}')
        flip_variables[(go_term, protein)] = model.NewBoolVar(f'flip_{go_term}_{protein}')
        proteins_added.add(protein)

    print(f"[INFO] Added {len(go_term_variables)} GO term variables and {len(flip_variables)} flip variables")
    print(f"[INFO] Added {len(proteins_added)} proteins")

    all_go_terms = sorted(all_go_terms) # sort for deterministic output

    # Add GO hierarchy constraints
    for child_go_term, protein in all_go_terms:
        for parent_go_term in sorted(go_hierarchy[child_go_term]):
            if (parent_go_term, protein) in all_go_terms:
                # child implies parent: child <= parent
                model.AddImplication(go_term_variables[(child_go_term, protein)], 
                                    go_term_variables[(parent_go_term, protein)])


    # Add heteromeric complex constraints
    for heteromeric_complex in sorted(participating_proteins.keys()):
        complex_vars = [go_term_variables[(heteromeric_complex, protein)] 
                        for protein in participating_proteins[heteromeric_complex]
                        if (heteromeric_complex, protein) in go_term_variables]
        if len(complex_vars) == 1:
            model.Add(complex_vars[0] == 0)
        elif complex_vars:
            no_members = model.NewBoolVar(f'complex_{heteromeric_complex}_has_zero_members')
            model.Add(sum(complex_vars) == 0).OnlyEnforceIf(no_members)
            model.Add(sum(complex_vars) >= 2).OnlyEnforceIf(no_members.Not())

    print(f"[INFO] Added {len(participating_proteins)} heteromeric complex constraints")

    # Add theta and eta constraints
    # Z3 formulation: theta = go_var XOR flip_var (for annotated terms)
    #                 eta = NOT(go_var) XOR flip_var (for non-annotated terms)
    # All theta and eta must be True (enforced via protein -> AND(thetas, etas) and genome -> AND(proteins))
    
    # Since genome and all proteins are forced True, we directly enforce:
    # For annotated terms: go_var XOR flip_var = True  =>  go_var + flip_var == 1
    # For non-annotated terms: NOT(go_var) XOR flip_var = True  =>  go_var == flip_var

    for go_term, protein in all_go_terms:
        if go_term in annotated_predictions[protein]:
            model.Add(go_term_variables[(go_term, protein)] + flip_variables[(go_term, protein)] == 1)
        else:
            model.Add(go_term_variables[(go_term, protein)] == flip_variables[(go_term, protein)])
    
    # for protein in sorted(annotated_predictions.keys()):
    #     # For annotated GO terms: XOR must be True
    #     for go_term in sorted(annotated_predictions[protein]):
    #         if (go_term, protein) in all_go_terms:
    #             # go_var XOR flip_var = True means exactly one is true
    #             model.Add(go_term_variables[(go_term, protein)] + flip_variables[(go_term, protein)] == 1)
        
    #     # For non-annotated GO terms: NOT(go_var) XOR flip_var = True
    #     for go_term, protein_pair in sorted(all_go_terms):
    #         if protein_pair == protein and go_term not in annotated_predictions[protein]:
    #             # NOT(go_var) XOR flip_var = True
    #             # This means: (1 - go_var) XOR flip_var = True
    #             # Which is: (1 - go_var) + flip_var = 1 (one of them is true)
    #             # Simplifies to: go_var = flip_var
    #             model.Add(go_term_variables[(go_term, protein)] == flip_variables[(go_term, protein)])

    print(f"[INFO] Added {len(all_go_terms)} flip conditions")

    # Add cost constraints - scale to integers with high precision
    # Use 1000000 scale factor to preserve precision down to 0.000001
    cost_scale = 1000000
    cost_terms = []
    for go_term, protein in all_go_terms:
        if (go_term, protein) in flip_variables:
            cost = int(flip_costs[protein][go_term] * cost_scale)
            # Ensure minimum cost of 1 to prevent "free" flips
            if cost == 0 and flip_costs[protein][go_term] > 0:
                cost = 1
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
        print(f"[INFO] SAT problem is unsatisfiable after {solve_time:.2f}s")
        print("[INFO] Returning original predictions")
        # Return original annotations as dict of bools
        result = defaultdict(dict)
        for protein, terms in annotated_predictions.items():
            for term in terms:
                result[protein][term] = True
        return result, 0


def get_participating_proteins(annotated_predictions: Dict[str, Set[str]], complex_terms: Set[str], predictions: Dict[str, Dict[str, float]], top_k = None, participating_threshold = None) -> Dict[str, List[str]]:
    """
    Get the participating proteins for each heteromeric complex.
    If top_k value is provided, consider only incoherent complexes and return the top-k proteins for each complex.
    If participating_threshold value is provided, return the proteins above the threshold for each complex.
    If neither are provided, return all proteins for each complex.
    """
    participating_proteins = defaultdict(list)
    if top_k is not None:
        for complex_go_term in complex_terms:
            # count how many proteins have the complex term
            proteins_with_complex = [p for p in annotated_predictions.keys() 
                                     if complex_go_term in annotated_predictions[p]]
            if len(proteins_with_complex) != 1:
                continue # only consider incoherent complexes
            
            complex_predictions_scores = [(p, predictions[p].get(complex_go_term, 0.0)) for p in predictions.keys()]
            sorted_proteins = sorted(complex_predictions_scores, key=lambda x: x[1], reverse=True)
            participating_proteins[complex_go_term] = [p for p, score in sorted_proteins[:top_k]]

            
    elif participating_threshold is not None:
        for complex_go_term in complex_terms:
            participating_proteins[complex_go_term] = [protein for protein in annotated_predictions.keys() 
                                                       if complex_go_term in predictions.get(protein, {}) 
                                                       and predictions[protein][complex_go_term] > participating_threshold]
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
    Get the final predictions from the adjusted predictions.
    """
    final_predictions = defaultdict(dict)
    for protein in predictions:
        for go_term in annotated_go_terms:
            if go_term not in predictions[protein]:
                final_score = 0
            else:
                final_score = predictions[protein][go_term]
            if go_term in adjusted_predictions[protein]:
                if adjusted_predictions[protein][go_term]:
                    final_score = max(final_score, fold_threshold + epsilon)
                else:
                    final_score = min(final_score, fold_threshold)
            if final_score > 0:
                final_predictions[protein][go_term] = final_score

    return final_predictions


def save_predictions(predictions: Dict[str, Dict[str, float]], file_path: str):
    """
    Save the adjusted predictions to a file.
    """
    with open(file_path, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        for protein, go_terms in predictions.items():
            go_terms_scores = [f"{go_term}|{score:.6f}" for go_term, score in go_terms.items()]

            writer.writerow([protein] + go_terms_scores)



def main(
    predictions_file: str,
    complexes_file: str,
    go_hierarchy_file: str,
    output_file: str,
    threshold: float = 0.5,
    optimized: bool = False,
    top_k: int = None,
    participating_threshold: float = None
):
    """
    Main function to run the complex coherence adjustment script.
    """
    
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
        participating_proteins = get_participating_proteins(annotated_predictions, heteromeric_complexes, predictions, top_k=top_k, participating_threshold=participating_threshold)

        print(f"[INFO] Number of Incoherent Complexes: {len(participating_proteins)}")

        if len(participating_proteins) == 0:
            print("[INFO] No participating proteins found, returning original predictions")
            if output_file != predictions_file:
                save_predictions(predictions, output_file)
            else:
                print("[INFO] Output file is the same as the predictions file, not saving")
            print("Done!")
            return 0, 0

        # compute after checking if there are participating proteins
        flip_costs, num_annotated_predictions = compute_flip_cost(predictions, threshold)
        
        go_protein_pairs = get_go_protein_pairs(annotated_predictions, go_hierarchy, participating_proteins)
        
        adjusted_predictions, total_flips = solve_sat_ortools_optimized(annotated_predictions, 
        go_protein_pairs,
        go_hierarchy, 
        flip_costs,
        participating_proteins)
    else:
        flip_costs, num_annotated_predictions = compute_flip_cost(predictions, threshold)

        adjusted_predictions, total_flips = solve_sat_ortools_hierarchy(annotated_predictions, 
        go_hierarchy, 
        flip_costs,
        heteromeric_complexes)

    print(f"[INFO] percentage of flipped predictions: {total_flips / num_annotated_predictions * 100:.2f}%")

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
    args = parser.parse_args()
    main(
        predictions_file=args.predictions,
        complexes_file=args.complexes,
        go_hierarchy_file=args.go_hierarchy,
        output_file=args.output,
        threshold=args.threshold,
        optimized=args.optimized,
        top_k=args.top_k,
        participating_threshold=args.participating_threshold
    )
