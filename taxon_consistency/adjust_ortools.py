"""
Adjust the predictions to be taxon consistent using min Flip SAT optimizer with OR-Tools.

Inputs:
- taxon_predictions.tsv: tab-separated file with protein name, GO term, and prediction score.
  Example: protein_name\tGO:term|score\tGO:term|score...
- go_taxon_constraints_updated.tsv: tab-separated file with GO term, taxon ID, and constraint type.
- go_hierarchy.tsv: tab-separated file with child GO term and parent GO term.

Outputs:
- taxon_adjusted_predictions.tsv: tab-separated file with protein name, GO term, and adjusted prediction score.

This version uses Google OR-Tools CP-SAT for 10-50x speedup over Z3.
"""

import csv
from collections import defaultdict
import itertools
from typing import Dict, List, Set, Tuple
from ortools.sat.python import cp_model
import time
import argparse


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

def load_constraints(file_path):
    """
    Load the constraints from the file.
    Expected format: GO_ID\tConstraint_Type\tTaxon_ID
    """
    in_taxon_constraints = defaultdict(list)
    never_in_taxon_constraints = defaultdict(list)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader, None)  # Skip header row
        for row in reader:
            if len(row) < 3:
                continue
            go_term = row[0]
            constraint_type = row[1]
            taxon_id = row[2]
            if constraint_type == 'only_in_taxon':
                in_taxon_constraints[go_term].append(taxon_id)
            elif constraint_type == 'never_in_taxon':
                never_in_taxon_constraints[go_term].append(taxon_id)
    return in_taxon_constraints, never_in_taxon_constraints

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


def load_taxon_hierarchy(file_path, ncbitaxon_hierarchy_file=None):
    """
    Load the taxon hierarchy from the file.
    Expected format:
    Term	Relationship	Parent/Disjoint_From_Term
    NCBITaxon_0000001	is_a	NCBITaxon_Union_0000000
    NCBITaxon_0000001	disjoint_from	NCBITaxon_Union_0000002
    NCBITaxon_Union_0000001	union_of	NCBITaxon_0000003
    """
    is_a_hierarchy = defaultdict(set)
    disjoint_from_hierarchy = defaultdict(set)
    union_of_hierarchy = defaultdict(set)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            taxon = row[0]
            relationship = row[1]
            parent_taxon = row[2]
            if relationship == 'is_a':
                is_a_hierarchy[taxon].add(parent_taxon)
            elif relationship == 'disjoint_from':
                disjoint_from_hierarchy[taxon].add(parent_taxon)
            elif relationship == 'union_of':
                union_of_hierarchy[taxon].add(parent_taxon)
    if ncbitaxon_hierarchy_file:
        with open(ncbitaxon_hierarchy_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                taxon = row[0]
                parent_taxon = row[1]
                is_a_hierarchy[taxon].add(parent_taxon)
    return is_a_hierarchy, disjoint_from_hierarchy, union_of_hierarchy


def compute_demotion_flip_cost(predictions: Dict[str, Dict[str, float]],
                       fold_threshold: float = 0.5) -> Dict[str, Dict[str, float]]:
    """
    Compute demotion flip cost for the predictions.
    """

    flip_cost_per_go_term = defaultdict(float)
    num_annotated_predictions = 0
    for protein in predictions:
        for go_term, score in predictions[protein].items():
            if score > fold_threshold:
                flip_cost_per_go_term[go_term] += score - fold_threshold
                num_annotated_predictions += 1

    print(f"Computed {len(flip_cost_per_go_term)} flip costs for {num_annotated_predictions} annotated predictions")
    return flip_cost_per_go_term, num_annotated_predictions


def get_annotated_predictions(predictions: Dict[str, Dict[str, float]], fold_threshold: float = 0.5) -> Dict[str, Set[str]]:
    """
    Get the annotated predictions from the predictions.
    """
    annotated_predictions = defaultdict(set)
    for protein in predictions:
        annotated_predictions[protein] = {
            go_term for go_term, score in predictions[protein].items() if score > fold_threshold
        }
    return annotated_predictions


def get_all_ancestors_of_list(taxon_hierarchy: Dict[str, Set[str]], taxon_list: List[str]) -> Set[str]:
    """
    Get the ancestors of a list of taxa from the taxon hierarchy.
    """
    all_ancestors = set()
    for taxon in taxon_list:
        all_ancestors.update(get_all_ancestors(taxon, taxon_hierarchy))
    return all_ancestors


def get_all_ancestors(taxon: str, hierarchy: Dict[str, Set[str]], 
                      cache: Dict[str, Set[str]] = None) -> Set[str]:
    """
    Get all ancestors of a taxon by traversing up the hierarchy.
    
    Includes the taxon itself in the result.
    Uses memoization for efficiency.
    """
    if cache is None:
        cache = {}
    
    if taxon in cache:
        return cache[taxon]
    
    ancestors = {taxon}  # Include the taxon itself
    parents = hierarchy.get(taxon, set())
    
    for parent in parents:
        ancestors.add(parent)
        # Recursively get ancestors of parent
        parent_ancestors = get_all_ancestors(parent, hierarchy, cache)
        ancestors.update(parent_ancestors)
    
    cache[taxon] = ancestors
    return ancestors


def normalize_taxon_id(taxon_id: str) -> str:
    """
    Normalize taxon ID to standard format.
    Handles both "83332" and "NCBITaxon_83332" formats.
    """
    if not taxon_id:
        return taxon_id
    
    # Remove NCBITaxon_ prefix if present
    if taxon_id.startswith('NCBITaxon_'):
        return taxon_id
    
    # Add NCBITaxon_ prefix if not present
    return f'NCBITaxon_{taxon_id}'

def solve_sat_ortools_hierarchy(annotated_predictions: Dict[str, Set[str]],
                in_taxon_constraints: Dict[str, List[str]],
                never_in_taxon_constraints: Dict[str, List[str]],
                go_hierarchy: Dict[str, Set[str]],
                flip_cost_per_go_term: Dict[str, float],
                taxon_is_a_hierarchy: Dict[str, Set[str]],
                taxon_disjoint_from_hierarchy: Dict[str, Set[str]],
                taxon_union_of_hierarchy: Dict[str, Set[str]],
                genome_taxon_id: str = None,
               ) -> Dict[str, Dict[str, bool]]:
    """
    Solve the SAT problem using OR-Tools CP-SAT.

    Args:
        annotated_predictions: Dict of sets of strings, where annotated_predictions[protein] is the set of go terms annotated to protein.
        in_taxon_constraints: Dict of lists of strings, where in_taxon_constraints[go_term] is the list of taxons that go_term is in.
        never_in_taxon_constraints: Dict of lists of strings, where never_in_taxon_constraints[go_term] is the list of taxons that go_term is never in.
        go_hierarchy: Dict of sets of strings, where go_hierarchy[go_term] is the set of parent go terms of go_term.
        flip_cost_per_go_term: Dict of floats, where flip_cost_per_go_term[go_term] is the cost of flipping the prediction for go_term.
        taxon_is_a_hierarchy: Dict of sets of strings, where taxon_is_a_hierarchy[taxon] is the set of parent taxons of taxon.
        taxon_disjoint_from_hierarchy: Dict of sets of strings, where taxon_disjoint_from_hierarchy[taxon] is the set of disjoint_from taxons of taxon.
        taxon_union_of_hierarchy: Dict of sets of strings, where taxon_union_of_hierarchy[taxon_union] is the set of union_of taxons of taxon_union (none of the children of taxon_union are union terms).
        genome_taxon_id: String, the taxon ID of the genome (e.g., "83332" or "NCBITaxon_83332"). If provided, constraints will be added to enforce that the genome belongs to this taxon.

    Returns:
        Dict of dicts of booleans, where adjusted_predictions[protein][go_term] is True if go_term is annotated to protein.
    """
    setup_start_time = time.time()
    model = cp_model.CpModel()

     # Add GO terms variables
    all_go_terms = set()
    go_term_variables = {}

    # Assumes annotations are propagated, so they should contain their parents.
    for protein in annotated_predictions:
        all_go_terms.update(annotated_predictions[protein])

    # Add GO term variables and flip variables
    flip_variables = {}
    for go_term in all_go_terms:
        go_term_variables[go_term] = model.NewBoolVar(f'go_{go_term}')
        flip_variables[go_term] = model.NewBoolVar(f'flip_{go_term}')
    print(f"[INFO] Added {len(go_term_variables)} GO term variables and {len(flip_variables)} flip variables")

    # Add GO hierarchy constraints
    for go_term in all_go_terms:
        for parent_go_term in go_hierarchy[go_term]:
            if parent_go_term in all_go_terms:
                model.AddImplication(go_term_variables[go_term], go_term_variables[parent_go_term])
    
    # get all taxons mapped by predictions
    in_taxon_taxons = set()
    never_in_taxon_taxons = set()
    for go_term in all_go_terms:
        in_taxon_taxons.update(in_taxon_constraints.get(go_term, []))
        never_in_taxon_taxons.update(never_in_taxon_constraints.get(go_term, []))
    common_taxons = in_taxon_taxons.intersection(never_in_taxon_taxons)
    all_taxons = in_taxon_taxons.union(never_in_taxon_taxons)
    print(f"[INFO] Common taxons count: {len(common_taxons)}")

    # Add a variable for each taxon 
    taxon_variables = {}
    for taxon in all_taxons:
        taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')
    print(f"[INFO] Added {len(taxon_variables)} taxon variables")

    # Add constraint for genome's taxon (if provided)
    if genome_taxon_id:
        normalized_genome_taxon = normalize_taxon_id(genome_taxon_id)
        if normalized_genome_taxon not in all_taxons:
            all_taxons.add(normalized_genome_taxon)
            taxon_variables[normalized_genome_taxon] = model.NewBoolVar(f'taxon_{normalized_genome_taxon}')
        model.Add(taxon_variables[normalized_genome_taxon] == 1)
        print(f"[INFO] Added constraint for genome taxon {genome_taxon_id}: {normalized_genome_taxon}")

    # performance optimization: only add taxon that are relevant
    relevant_taxon_is_a_hierarchy = get_all_ancestors_of_list(taxon_is_a_hierarchy, all_taxons)
    # add relevant taxon is_a hierarchy to taxon_variables that are not already in taxon_variables
    for taxon in relevant_taxon_is_a_hierarchy:
        if taxon not in all_taxons:
            all_taxons.add(taxon)
            taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')

    relevant_taxon_disjoint_from_hierarchy = get_all_ancestors_of_list(taxon_disjoint_from_hierarchy, all_taxons)
    # add relevant taxon disjoint_from hierarchy to taxon_variables that are not already in taxon_variables
    for taxon in relevant_taxon_disjoint_from_hierarchy:
        if taxon not in all_taxons:
            all_taxons.add(taxon)
            taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')


    # Add taxon is_a hierarchy constraints
    count_is_a_constraints = 0
    for taxon in all_taxons:
        for parent_taxon in taxon_is_a_hierarchy[taxon]:
            if parent_taxon in all_taxons:

                model.AddImplication(taxon_variables[taxon], taxon_variables[parent_taxon])
                count_is_a_constraints += 1

    # Add taxon disjoint_from hierarchy constraints
    count_disjoint_from_constraints = 0
    for taxon in all_taxons:
        for disjoint_from_taxon in taxon_disjoint_from_hierarchy[taxon]:
            if disjoint_from_taxon in all_taxons:

                model.AddBoolOr([taxon_variables[taxon].Not(), taxon_variables[disjoint_from_taxon].Not()])
                count_disjoint_from_constraints += 1

    print(f"[INFO] Added {count_is_a_constraints} is_a constraints")
    print(f"[INFO] Added {count_disjoint_from_constraints} disjoint_from constraints")


    # Add go and taxon variables constraints
    taxons_added_with_go_term_constraints = set()    
    for go_term in all_go_terms:
        for taxon in in_taxon_constraints[go_term]:
            if taxon in taxon_variables:
                model.AddImplication(go_term_variables[go_term], taxon_variables[taxon])
                taxons_added_with_go_term_constraints.add(taxon)
    for go_term in all_go_terms:
        for taxon in never_in_taxon_constraints[go_term]:
            if taxon in taxon_variables:
                model.AddImplication(go_term_variables[go_term], taxon_variables[taxon].Not())
                taxons_added_with_go_term_constraints.add(taxon)

    # Add taxon union constraints (IMPORTANT: This must be added after the taxon variables are added)
    disjoint_children_constraints = 0
    for taxon in all_taxons:
        if taxon.startswith('NCBITaxon_Union_'):
            children_taxons = [child_taxon for child_taxon in all_taxons if taxon in taxon_is_a_hierarchy[child_taxon]]
            # add union equivalent to OR of children_taxons
            if children_taxons:
                children_vars = [taxon_variables[child_taxon] for child_taxon in children_taxons]
                # Union taxon is True if and only if at least one child is True
                # taxon == OR(children) means: taxon -> OR(children) AND OR(children) -> taxon
                # If taxon is True, at least one child must be True
                model.AddBoolOr(children_vars).OnlyEnforceIf(taxon_variables[taxon])
                # If any child is True, taxon must be True
                for child_var in children_vars:
                    model.AddImplication(child_var, taxon_variables[taxon])

            # This is currently only adding disjoint constraints for union children. TODO: check if this is sufficient.
            if taxon in taxon_union_of_hierarchy:
                taxon_union_of_children = taxon_union_of_hierarchy[taxon]
                taxon_union_of_children = taxon_union_of_children.intersection(all_taxons)
                if len(taxon_union_of_children) > 1:
                    for child_taxon1, child_taxon2 in itertools.combinations(taxon_union_of_children, 2):
                        model.AddBoolOr([taxon_variables[child_taxon1].Not(), taxon_variables[child_taxon2].Not()])
                        disjoint_children_constraints += 1
        else:
            # if taxon in taxons_added_with_go_term_constraints:
            children_taxons = [child_taxon for child_taxon in all_taxons if taxon in taxon_is_a_hierarchy[child_taxon]]
            # filter out union taxons. TODO: I don't think this is complete. Might want to replace union taxons with their children, but this is not trivial as a child might also be in children_taxons already. Can't have not(And(child_1, child_1)).
            children_taxons = [child_taxon for child_taxon in children_taxons if not child_taxon.startswith('NCBITaxon_Union_')]

            if len(children_taxons) > 1:
                for child_taxon1, child_taxon2 in itertools.combinations(children_taxons, 2):
                    model.AddBoolOr([taxon_variables[child_taxon1].Not(), taxon_variables[child_taxon2].Not()])
                    disjoint_children_constraints += 1

    print(f"[INFO] Added {disjoint_children_constraints} disjoint children constraints")

    

    # Add XOR constraints for each GO term (enforce once per GO term, not per protein)
    # XOR constraint: go_var XOR flip_var = True means exactly one is true: go_var + flip_var == 1
    for go_term in all_go_terms:
        model.Add(go_term_variables[go_term] + flip_variables[go_term] == 1)

    # No need to add protein variables and constraints since we are enforcing the XOR constraints globally.
    # Similarly, no need to add genome variable and constraints since we are enforcing the XOR constraints globally.

    # Add cost constraints - scale to integers with high precision
    # Use 1000000 scale factor to preserve precision down to 0.000001
    cost_scale = 1000000

    def compute_cost(go_term):
        cost = int(flip_cost_per_go_term[go_term] * cost_scale)
        return cost if cost > 0 else (1 if flip_cost_per_go_term[go_term] > 0 else 0)

    cost_terms = [
        flip_variables[go_term] * compute_cost(go_term)
        for go_term in flip_variables
    ]
    
    if cost_terms:
        model.Minimize(sum(cost_terms))
    else:
        # If no costs, just minimize number of flips
        # Convert boolean variables to integers (1 if True, 0 if False) for summing
        flip_sum = sum([flip_variables[go_term] for go_term in flip_variables])
        model.Minimize(flip_sum)

    solve_start_time = time.time()
    print(f"[INFO] Setup time: {solve_start_time - setup_start_time:.2f}s")
    print("[INFO] Solving SAT problem...")
    
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    solve_time = time.time() - solve_start_time


    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adjusted_predictions = defaultdict(dict)
        print(f"[INFO] SAT problem solved successfully in {solve_time:.2f}s")
        
        # Pre-compute all flip results once (much faster than computing in nested loop)
        flipped_go_terms = {
            go_term for go_term in flip_variables 
            if solver.BooleanValue(flip_variables[go_term])
        }
        
        # Build adjusted predictions using dictionary comprehension (faster than repeated dict assignments)
        total_flips = 0
        for protein, go_terms in annotated_predictions.items():
            # Count flips for this protein and build the dict in one pass
            protein_flips = flipped_go_terms & go_terms
            total_flips += len(protein_flips)
            
            # Dictionary comprehension: False only if in flipped set, True otherwise
            adjusted_predictions[protein] = {
                go_term: (go_term not in flipped_go_terms)
                for go_term in go_terms
            }
        
        print(f"[INFO]   Total flips: {total_flips}")

        
        opt_cost = solver.ObjectiveValue() / cost_scale if cost_terms else 0.0
        print(f"[INFO]   Optimal cost: {opt_cost}")

        return adjusted_predictions, total_flips
    else:
        print(f"[INFO] SAT problem is unsatisfiable after {solve_time:.2f}s")
        print(f"[INFO] Returning original predictions")    
        # Return original annotations as dict of bools
        result = defaultdict(dict)
        for protein, terms in annotated_predictions.items():
            for term in terms:
                result[protein][term] = True
        return result, 0


def solve_sat_ortools_taxon_assignment(annotated_predictions: Dict[str, Set[str]],
                in_taxon_constraints: Dict[str, List[str]],
                never_in_taxon_constraints: Dict[str, List[str]],
                go_hierarchy: Dict[str, Set[str]],
                flip_cost_per_go_term: Dict[str, float],
                taxon_is_a_hierarchy: Dict[str, Set[str]],
                taxon_disjoint_from_hierarchy: Dict[str, Set[str]],
                taxon_union_of_hierarchy: Dict[str, Set[str]],
                genome_taxon_id: str = None,
               ) -> Tuple[Dict[str, Dict[str, bool]], List[str], int]:
    """
    Solve the SAT problem using OR-Tools CP-SAT and obtain the predicted taxon assignment.

    Args:
        annotated_predictions: Dict of sets of strings, where annotated_predictions[protein] is the set of go terms annotated to protein.
        in_taxon_constraints: Dict of lists of strings, where in_taxon_constraints[go_term] is the list of taxons that go_term is in.
        never_in_taxon_constraints: Dict of lists of strings, where never_in_taxon_constraints[go_term] is the list of taxons that go_term is never in.
        go_hierarchy: Dict of sets of strings, where go_hierarchy[go_term] is the set of parent go terms of go_term.
        flip_cost_per_go_term: Dict of floats, where flip_cost_per_go_term[go_term] is the cost of flipping the prediction for go_term.
        taxon_is_a_hierarchy: Dict of sets of strings, where taxon_is_a_hierarchy[taxon] is the set of parent taxons of taxon.
        taxon_disjoint_from_hierarchy: Dict of sets of strings, where taxon_disjoint_from_hierarchy[taxon] is the set of disjoint_from taxons of taxon.
        taxon_union_of_hierarchy: Dict of sets of strings, where taxon_union_of_hierarchy[taxon_union] is the set of union_of taxons of taxon_union (none of the children of taxon_union are union terms).
        genome_taxon_id: String, the taxon ID of the genome (e.g., "83332" or "NCBITaxon_83332"). If provided, constraints will be added to enforce that the genome belongs to this taxon.

    Returns:
        Dict of dicts of booleans, where adjusted_predictions[protein][go_term] is True if go_term is annotated to protein.
        List of strings, the predicted taxon assignment.
        int, the total number of flips.
    """
    setup_start_time = time.time()
    model = cp_model.CpModel()

    # Add GO terms variables
    all_go_terms = set()
    go_term_variables = {}

    # Assumes annotations are propagated, so they should contain their parents.
    for protein in annotated_predictions:
        all_go_terms.update(annotated_predictions[protein])

    # Add GO term variables and flip variables
    flip_variables = {}
    for go_term in all_go_terms:
        go_term_variables[go_term] = model.NewBoolVar(f'go_{go_term}')
        flip_variables[go_term] = model.NewBoolVar(f'flip_{go_term}')
    print(f"[INFO] Added {len(go_term_variables)} GO term variables and {len(flip_variables)} flip variables")

    # Add GO hierarchy constraints
    for go_term in all_go_terms:
        for parent_go_term in go_hierarchy[go_term]:
            if parent_go_term in all_go_terms:
                model.AddImplication(go_term_variables[go_term], go_term_variables[parent_go_term])
    
    # get all taxons mapped by predictions
    in_taxon_taxons = set()
    never_in_taxon_taxons = set()
    for go_term in all_go_terms:
        in_taxon_taxons.update(in_taxon_constraints.get(go_term, []))
        never_in_taxon_taxons.update(never_in_taxon_constraints.get(go_term, []))
    common_taxons = in_taxon_taxons.intersection(never_in_taxon_taxons)
    all_taxons = in_taxon_taxons.union(never_in_taxon_taxons)
    print(f"[INFO] Common taxons count: {len(common_taxons)}")

    # Add a variable for each taxon 
    taxon_variables = {}
    for taxon in all_taxons:
        taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')
    print(f"[INFO] Added {len(taxon_variables)} taxon variables")

    # Add constraint for genome's taxon (if provided)
    if genome_taxon_id:
        normalized_genome_taxon = normalize_taxon_id(genome_taxon_id)
        if normalized_genome_taxon not in all_taxons:
            all_taxons.add(normalized_genome_taxon)
            taxon_variables[normalized_genome_taxon] = model.NewBoolVar(f'taxon_{normalized_genome_taxon}')
        model.Add(taxon_variables[normalized_genome_taxon] == 1)
        print(f"[INFO] Added constraint for genome taxon {genome_taxon_id}: {normalized_genome_taxon}")

    # performance optimization: only add taxon that are relevant
    relevant_taxon_is_a_hierarchy = get_all_ancestors_of_list(taxon_is_a_hierarchy, all_taxons)
    # add relevant taxon is_a hierarchy to taxon_variables that are not already in taxon_variables
    for taxon in relevant_taxon_is_a_hierarchy:
        if taxon not in all_taxons:
            all_taxons.add(taxon)
            taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')

    relevant_taxon_disjoint_from_hierarchy = get_all_ancestors_of_list(taxon_disjoint_from_hierarchy, all_taxons)
    # add relevant taxon disjoint_from hierarchy to taxon_variables that are not already in taxon_variables
    for taxon in relevant_taxon_disjoint_from_hierarchy:
        if taxon not in all_taxons:
            all_taxons.add(taxon)
            taxon_variables[taxon] = model.NewBoolVar(f'taxon_{taxon}')

    # Add taxon is_a hierarchy constraints
    count_is_a_constraints = 0
    for taxon in all_taxons:
        for parent_taxon in taxon_is_a_hierarchy[taxon]:
            if parent_taxon in all_taxons:
                model.AddImplication(taxon_variables[taxon], taxon_variables[parent_taxon])
                count_is_a_constraints += 1

    # Add taxon disjoint_from hierarchy constraints
    count_disjoint_from_constraints = 0
    for taxon in all_taxons:
        for disjoint_from_taxon in taxon_disjoint_from_hierarchy[taxon]:
            if disjoint_from_taxon in all_taxons:
                model.AddBoolOr([taxon_variables[taxon].Not(), taxon_variables[disjoint_from_taxon].Not()])
                count_disjoint_from_constraints += 1
    print(f"[INFO] Added {count_is_a_constraints} is_a constraints")
    print(f"[INFO] Added {count_disjoint_from_constraints} disjoint_from constraints")

    # Add go and taxon variables constraints
    taxons_added_with_go_term_constraints = set()    
    for go_term in all_go_terms:
        for taxon in in_taxon_constraints[go_term]:
            if taxon in taxon_variables:
                model.AddImplication(go_term_variables[go_term], taxon_variables[taxon])
                taxons_added_with_go_term_constraints.add(taxon)
    for go_term in all_go_terms:
        for taxon in never_in_taxon_constraints[go_term]:
            if taxon in taxon_variables:
                model.AddImplication(go_term_variables[go_term], taxon_variables[taxon].Not())
                taxons_added_with_go_term_constraints.add(taxon)

    # Add taxon union constraints (IMPORTANT: This must be added after the taxon variables are added)
    disjoint_children_constraints = 0
    for taxon in all_taxons:
        if taxon.startswith('NCBITaxon_Union_'):
            children_taxons = [child_taxon for child_taxon in all_taxons if taxon in taxon_is_a_hierarchy[child_taxon]]
            # add union equivalent to OR of children_taxons
            if children_taxons:
                children_vars = [taxon_variables[child_taxon] for child_taxon in children_taxons]
                # Union taxon is True if and only if at least one child is True
                # taxon == OR(children) means: taxon -> OR(children) AND OR(children) -> taxon
                # If taxon is True, at least one child must be True
                model.AddBoolOr(children_vars).OnlyEnforceIf(taxon_variables[taxon])
                # If any child is True, taxon must be True
                for child_var in children_vars:
                    model.AddImplication(child_var, taxon_variables[taxon])
            # This is currently only adding disjoint constraints for union children. TODO: check if this is sufficient.
            if taxon in taxon_union_of_hierarchy:
                taxon_union_of_children = taxon_union_of_hierarchy[taxon]
                taxon_union_of_children = taxon_union_of_children.intersection(all_taxons)
                if len(taxon_union_of_children) > 1:
                    for child_taxon1, child_taxon2 in itertools.combinations(taxon_union_of_children, 2):
                        model.AddBoolOr([taxon_variables[child_taxon1].Not(), taxon_variables[child_taxon2].Not()])
                        disjoint_children_constraints += 1
        else:
            # if taxon in taxons_added_with_go_term_constraints:
            children_taxons = [child_taxon for child_taxon in all_taxons if taxon in taxon_is_a_hierarchy[child_taxon]]
            # filter out union taxons. TODO: I don't think this is complete. Might want to replace union taxons with their children, but this is not trivial as a child might also be in children_taxons already. Can't have not(And(child_1, child_1)).
            children_taxons = [child_taxon for child_taxon in children_taxons if not child_taxon.startswith('NCBITaxon_Union_')]
            if len(children_taxons) > 1:
                for child_taxon1, child_taxon2 in itertools.combinations(children_taxons, 2):
                    model.AddBoolOr([taxon_variables[child_taxon1].Not(), taxon_variables[child_taxon2].Not()])
                    disjoint_children_constraints += 1

    print(f"[INFO] Added {disjoint_children_constraints} disjoint children constraints")

    for go_term in all_go_terms:
        model.Add(go_term_variables[go_term] + flip_variables[go_term] == 1)
    cost_scale = 1000000

    def compute_cost(go_term):
        cost = int(flip_cost_per_go_term[go_term] * cost_scale)
        return cost if cost > 0 else (1 if flip_cost_per_go_term[go_term] > 0 else 0)

    cost_terms = [
        flip_variables[go_term] * compute_cost(go_term)
        for go_term in flip_variables
    ]
    if cost_terms:
        model.Minimize(sum(cost_terms))
    else:
        # If no costs, just minimize number of flips
        # Convert boolean variables to integers (1 if True, 0 if False) for summing
        flip_sum = sum([flip_variables[go_term] for go_term in flip_variables])
        model.Minimize(flip_sum)

    solve_start_time = time.time()
    print(f"[INFO] Setup time: {solve_start_time - setup_start_time:.2f}s")
    print("[INFO] Solving SAT problem...")
    solver = cp_model.CpSolver()
    status = solver.Solve(model)
    solve_time = time.time() - solve_start_time

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adjusted_predictions = defaultdict(dict)
        print(f"[INFO] SAT problem solved successfully in {solve_time:.2f}s")
        # Pre-compute all flip results once (much faster than computing in nested loop)
        flipped_go_terms = {
            go_term for go_term in flip_variables 
            if solver.BooleanValue(flip_variables[go_term])
        }
        
        # Build adjusted predictions using dictionary comprehension (faster than repeated dict assignments)
        total_flips = 0
        for protein, go_terms in annotated_predictions.items():
            # Count flips for this protein and build the dict in one pass
            protein_flips = flipped_go_terms & go_terms
            total_flips += len(protein_flips)
            # Dictionary comprehension: False only if in flipped set, True otherwise
            adjusted_predictions[protein] = {
                go_term: (go_term not in flipped_go_terms)
                for go_term in go_terms
            }
        print(f"[INFO]   Total flips: {total_flips}")

        predicted_taxons = [taxon for taxon in taxon_variables if solver.BooleanValue(taxon_variables[taxon])]
        print(f"[INFO]   taxon assignments: {len(predicted_taxons)}")
        opt_cost = solver.ObjectiveValue() / cost_scale if cost_terms else 0.0
        print(f"[INFO]   Optimal cost: {opt_cost}")

        return adjusted_predictions, predicted_taxons, total_flips
    else:
        print(f"[INFO] SAT problem is unsatisfiable after {solve_time:.2f}s")
        print(f"[INFO] Returning original predictions")    
        # Return original annotations as dict of bools
        result = defaultdict(dict)
        for protein, terms in annotated_predictions.items():
            for term in terms:
                result[protein][term] = True
        return result, [], 0


def get_final_predictions(adjusted_predictions: Dict[str, Dict[str, bool]], predictions: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """
    Get the final predictions from the adjusted predictions.
    Remove the demoted predictions entirely.
    """
    # Use dictionary comprehension with .get() for fast filtering
    # Filter out the demoted predictions entirely.
    final_predictions = {
        protein: {
            go_term: score 
            for go_term, score in go_terms.items()
            if adjusted_predictions.get(protein, {}).get(go_term, True)
        }
        for protein, go_terms in predictions.items()
    }
    
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


def adjust_per_taxon(
    predictions_file: str,
    in_taxon_constraints: Dict[str, List[str]],
    never_in_taxon_constraints: Dict[str, List[str]],
    go_hierarchy: Dict[str, Set[str]],
    taxon_is_a_hierarchy: Dict[str, Set[str]],
    taxon_disjoint_from_hierarchy: Dict[str, Set[str]],
    taxon_union_of_hierarchy: Dict[str, Set[str]],
    output_file: str,
    threshold: float = 0.5,
    taxon_id: str = None
):
    predictions = load_predictions(predictions_file)
    print(f"Loaded {len(predictions)} proteins for taxon {taxon_id}")

    annotated_predictions = get_annotated_predictions(predictions, threshold)
    flip_cost_per_go_term, num_annotated_predictions = compute_demotion_flip_cost(predictions, threshold)

    print("Setting up SAT problem...")
    adjusted_predictions, total_flips = solve_sat_ortools_hierarchy(annotated_predictions, 
        in_taxon_constraints, 
        never_in_taxon_constraints, 
        go_hierarchy, 
        flip_cost_per_go_term,
        taxon_is_a_hierarchy,
        taxon_disjoint_from_hierarchy,
        taxon_union_of_hierarchy,
        taxon_id,
        )

    print(f"[INFO] percentage of flipped predictions: {total_flips / num_annotated_predictions * 100:.2f}%")

    
    final_predictions = get_final_predictions(adjusted_predictions, predictions)
    
    print(f"\nSaving final predictions to {output_file}...")
    save_predictions(final_predictions, output_file)
    print("Done!")

    return total_flips, num_annotated_predictions

def main(
    predictions_file: str,
    constraints_file: str,
    go_hierarchy_file: str,
    taxon_hierarchy_file: str,
    ncbitaxon_hierarchy_file: str,
    output_file: str,
    threshold: float = 0.5,
    taxon_id: str = None
):
    """
    Main function to run the taxon adjustment script.
    """
    print("Loading data...")
    predictions = load_predictions(predictions_file)
    in_taxon_constraints, never_in_taxon_constraints = load_constraints(constraints_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file) if go_hierarchy_file else defaultdict(set)
    taxon_is_a_hierarchy, taxon_disjoint_from_hierarchy, taxon_union_of_hierarchy = load_taxon_hierarchy(taxon_hierarchy_file, ncbitaxon_hierarchy_file) if taxon_hierarchy_file else (defaultdict(set), defaultdict(set), defaultdict(set))

    print(f"Loaded {len(predictions)} proteins")
    print(f"Loaded {len(in_taxon_constraints)} 'only_in_taxon' constraints")
    print(f"Loaded {len(never_in_taxon_constraints)} 'never_in_taxon' constraints")
    print(f"Loaded {len(go_hierarchy)} GO hierarchy relationships")
    print(f"Loaded {len(taxon_is_a_hierarchy)} taxon is_a hierarchy relationships")
    print(f"Loaded {len(taxon_disjoint_from_hierarchy)} taxon disjoint_from hierarchy relationships")
    print(f"Loaded {len(taxon_union_of_hierarchy)} taxon union_of hierarchy relationships")


    annotated_predictions = get_annotated_predictions(predictions, threshold)

    flip_cost_per_go_term, num_annotated_predictions = compute_demotion_flip_cost(predictions, threshold)

    print("Setting up SAT problem...")
    adjusted_predictions, total_flips = solve_sat_ortools_hierarchy(annotated_predictions, 
        in_taxon_constraints, 
        never_in_taxon_constraints, 
        go_hierarchy, 
        flip_cost_per_go_term,
        taxon_is_a_hierarchy,
        taxon_disjoint_from_hierarchy,
        taxon_union_of_hierarchy,
        taxon_id,
        )

    print(f"[INFO] percentage of flipped predictions: {total_flips / num_annotated_predictions * 100:.2f}%")

    final_predictions = get_final_predictions(adjusted_predictions, predictions)
    
    print(f"\nSaving final predictions to {output_file}...")
    save_predictions(final_predictions, output_file)
    print("Done!")

    return total_flips, num_annotated_predictions


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Adjust predictions to be taxon consistent using min Flip SAT optimizer with OR-Tools.')
    parser.add_argument('--predictions', required=True, help='Path to taxon_predictions.tsv')
    parser.add_argument('--constraints', required=True, help='Path to go_taxon_constraints_updated.tsv')
    parser.add_argument('--go-hierarchy', required=True, help='Path to go_hierarchy.tsv')
    parser.add_argument('--taxon-hierarchy', required=True, help='Path to taxon_hierarchy.tsv')
    parser.add_argument('--output', default=None, help='Path to output taxon_adjusted_predictions.tsv')
    parser.add_argument('--threshold', type=float, default=0.5, help='Prediction threshold (default: 0.5)')
    parser.add_argument('--ncbitaxon-hierarchy', required=True, help='Path to ncbitaxon_hierarchy.tsv')
    parser.add_argument('--taxon-id', required=False, default=None, help='Taxon ID of the genome (e.g., 83332). If not provided, will try to extract from filename.')
    
    args = parser.parse_args()
    if args.output is None:
        # if predictions path has */predictions/predictions_*, then output path should be */optimized/optimized_*
        if 'predictions/' in args.predictions:
            args.output = args.predictions.replace('predictions/', 'optimized/')
            args.output = args.output.replace('predictions_', 'optimized_')
        else:
            args.output = args.predictions.replace('.tsv', '_adjusted.tsv')
    main(
        predictions_file=args.predictions,
        constraints_file=args.constraints,
        go_hierarchy_file=args.go_hierarchy,
        taxon_hierarchy_file=args.taxon_hierarchy,
        output_file=args.output,
        threshold=args.threshold,
        ncbitaxon_hierarchy_file=args.ncbitaxon_hierarchy,
        taxon_id=args.taxon_id
    )

