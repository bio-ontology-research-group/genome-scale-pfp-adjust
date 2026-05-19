"""
This script finds taxon inconsistencies in test set protein annotations.

Combines patterns from:
- collect_preds.py - loading test proteins by taxon from TSV files
- check_swissprot_annots_taxon_consistency.py - checking taxon consistency

Detects GO term annotations that violate taxon constraints:
- only_in_taxon: GO term should only appear in specific taxa
- never_in_taxon: GO term should never appear in specific taxa
"""

import csv
import os
import sys
import glob
import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
from tqdm import tqdm

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Default paths
DEFAULT_ANNOTATIONS_DIR = "${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic"
DEFAULT_CONSTRAINTS = os.path.join(BASE_DIR, 'data', 'go_taxon_constraints_extracted_obo.tsv')
DEFAULT_GO_HIERARCHY = os.path.join(BASE_DIR, 'data', 'go_hierarchy.tsv')
DEFAULT_TAXON_HIERARCHY = os.path.join(BASE_DIR, 'data', 'taxon_hierarchy.tsv')
DEFAULT_NCBITAXON_HIERARCHY = os.path.join(BASE_DIR, 'data', 'ncbitaxon_hierarchy.tsv')
DEFAULT_GO_FILE = os.path.join(BASE_DIR, 'data', 'go-basic.obo')


def load_test_proteins(test_proteins_file: str) -> Dict[str, List[str]]:
    """
    Load the test proteins from the test proteins file with matching taxon IDs.

    Example: (proteins_by_date_23-MAY-2024.tsv file)
    protein_id	protein_name	integration_date	taxon_id	taxon_id_matched
    10H_STRNX	Strychnine-10-hydroxylase {ECO:0000303|PubMed:35794473}	27-NOV-2024	NA	False
    34K1_AEDAL	Salivary protein SG34 {ECO:0000305}	02-OCT-2024	7160	True
    ...

    Returns:
        Dict[str, List[str]]: Dictionary mapping taxon IDs to list of test protein IDs.
    """
    test_proteins = defaultdict(list)
    with open(test_proteins_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            protein_id = parts[0]
            if protein_id == 'protein_id':
                continue
            taxon_id = parts[3]
            if taxon_id != 'NA':
                test_proteins[taxon_id].append(protein_id)
    return test_proteins


def load_annotations(
    annotations_dir: str,
    test_proteins: Dict[str, List[str]],
    subontology: Optional[str] = None,
    go=None
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Load GO annotations for test proteins from per-taxon annotation files.

    Annotation files have no header and format: protein_id\tGO:term\tGO:term\t...

    Args:
        annotations_dir: Directory containing annots_taxon_{taxon_id}.tsv files
        test_proteins: Dict mapping taxon IDs to list of protein IDs
        subontology: Optional filter for GO subontology (cc, bp, mf)
        go: Ontology object for namespace filtering (required if subontology is set)

    Returns:
        Dict mapping organism_id -> Dict mapping protein_id -> set of GO terms
    """
    namespace_map = {
        'cc': 'cellular_component',
        'bp': 'biological_process',
        'mf': 'molecular_function'
    }

    annotations = defaultdict(dict)

    for taxon_id, protein_ids in tqdm(test_proteins.items(), desc="Loading annotations"):
        # Find annotation file for this taxon
        annotation_files = glob.glob(os.path.join(annotations_dir, f'annots_taxon_{taxon_id}.tsv'))
        if len(annotation_files) == 0:
            continue
        annotation_file = annotation_files[0]

        protein_id_set = set(protein_ids)
        with open(annotation_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                protein_id = parts[0]
                if protein_id not in protein_id_set:
                    continue

                go_terms = set(parts[1:])

                # Filter by subontology if specified
                if subontology and go:
                    target_namespace = namespace_map.get(subontology)
                    if target_namespace:
                        filtered_terms = set()
                        for term in go_terms:
                            if go.has_term(term):
                                if go.get_namespace(term) == target_namespace:
                                    filtered_terms.add(term)
                        if filtered_terms:
                            annotations[taxon_id][protein_id] = filtered_terms

    return annotations


def load_constraints(file_path: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
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


def load_go_hierarchy(file_path: str) -> Dict[str, Set[str]]:
    """
    Load the GO hierarchy from the file.
    Format: child_GO_term\tparent_GO_term
    """
    hierarchy = defaultdict(set)
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue
            child_go_term = row[0]
            parent_go_term = row[1]
            hierarchy[child_go_term].add(parent_go_term)
    return hierarchy


def load_taxon_hierarchy(
    file_path: str,
    ncbitaxon_hierarchy_file: Optional[str] = None
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
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
            if len(row) < 3:
                continue
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
                if len(row) < 2:
                    continue
                taxon = row[0]
                parent_taxon = row[1]
                is_a_hierarchy[taxon].add(parent_taxon)
    return is_a_hierarchy, disjoint_from_hierarchy, union_of_hierarchy


def get_all_ancestors(
    taxon: str,
    hierarchy: Dict[str, Set[str]],
    cache: Optional[Dict[str, Set[str]]] = None,
    visited: Optional[Set[str]] = None
) -> Set[str]:
    """
    Get all ancestors of a taxon by traversing up the hierarchy.

    Includes the taxon itself in the result.
    Uses memoization for efficiency and cycle detection to prevent infinite recursion.
    """
    if cache is None:
        cache = {}
    if visited is None:
        visited = set()

    # If already computed, return cached result
    if taxon in cache:
        return cache[taxon]

    # If currently being processed (cycle detected), return taxon itself to break cycle
    if taxon in visited:
        return {taxon}

    # Mark as being processed
    visited.add(taxon)

    ancestors = {taxon}  # Include the taxon itself
    parents = hierarchy.get(taxon, set())

    for parent in parents:
        ancestors.add(parent)
        # Recursively get ancestors of parent
        parent_ancestors = get_all_ancestors(parent, hierarchy, cache, visited)
        ancestors.update(parent_ancestors)

    # Mark as processed and cache the complete result
    visited.remove(taxon)
    cache[taxon] = ancestors
    return ancestors


def normalize_taxon_id(taxon_id: str) -> str:
    """
    Normalize taxon ID to standard format.
    Handles both "83332" and "NCBITaxon_83332" formats.
    """
    if not taxon_id:
        return taxon_id

    # If already has NCBITaxon_ prefix, return as is
    if taxon_id.startswith('NCBITaxon_'):
        return taxon_id

    # Add NCBITaxon_ prefix if not present
    return f'NCBITaxon_{taxon_id}'


def get_all_go_ancestors(
    go_term: str,
    go_hierarchy: Dict[str, Set[str]],
    cache: Optional[Dict[str, Set[str]]] = None,
    visited: Optional[Set[str]] = None
) -> Set[str]:
    """
    Get all ancestors of a GO term by traversing up the hierarchy.

    Includes the GO term itself in the result.
    Uses memoization for efficiency and cycle detection to prevent infinite recursion.
    """
    if cache is None:
        cache = {}
    if visited is None:
        visited = set()

    # If already computed, return cached result
    if go_term in cache:
        return cache[go_term]

    # If currently being processed (cycle detected), return GO term itself to break cycle
    if go_term in visited:
        return {go_term}

    # Mark as being processed
    visited.add(go_term)

    ancestors = {go_term}  # Include the GO term itself
    parents = go_hierarchy.get(go_term, set())

    for parent in parents:
        ancestors.add(parent)
        # Recursively get ancestors of parent
        parent_ancestors = get_all_go_ancestors(parent, go_hierarchy, cache, visited)
        ancestors.update(parent_ancestors)

    # Mark as processed and cache the complete result
    visited.remove(go_term)
    cache[go_term] = ancestors
    return ancestors


def check_consistency(
    annotations_by_organism: Dict[str, Dict[str, Set[str]]],
    in_taxon_constraints: Dict[str, List[str]],
    never_in_taxon_constraints: Dict[str, List[str]],
    taxon_is_a_hierarchy: Optional[Dict[str, Set[str]]] = None,
    go_hierarchy: Optional[Dict[str, Set[str]]] = None
) -> List[Tuple[str, str, str, str, str]]:
    """
    Check taxon consistency for all protein-GO term pairs.

    Args:
        annotations_by_organism: Dict mapping organism ID -> Dict mapping protein ID -> set of GO terms
        in_taxon_constraints: Dict mapping GO term -> list of allowed taxon IDs
        never_in_taxon_constraints: Dict mapping GO term -> list of forbidden taxon IDs
        taxon_is_a_hierarchy: Dict mapping taxon -> set of parent taxa (optional)
        go_hierarchy: Dict mapping GO term -> set of parent GO terms (optional)

    Returns:
        List of tuples: (protein_id, GO_term, organism_id, constraint_type,
                        conflicting_taxon_constraint)
    """
    go_ancestor_cache = {}
    inconsistencies = []
    inconsistent_organisms = set()

    for organism_id in tqdm(annotations_by_organism, desc="Checking consistency for all proteins"):
        normalized_taxon = normalize_taxon_id(organism_id)

        # Build ancestor cache for taxon hierarchy
        taxon_ancestor_cache = {}
        if taxon_is_a_hierarchy:
            organism_ancestors = get_all_ancestors(normalized_taxon, taxon_is_a_hierarchy, taxon_ancestor_cache)
        else:
            organism_ancestors = {normalized_taxon}

        for protein_id in annotations_by_organism[organism_id]:
            for go_term in annotations_by_organism[organism_id][protein_id]:
                # Check constraints for this GO term and all its ancestors
                go_terms_to_check = {go_term}
                if go_hierarchy:
                    go_terms_to_check.update(get_all_go_ancestors(go_term, go_hierarchy, go_ancestor_cache))

                # Check never_in_taxon constraints
                never_in_violation_found = False
                for check_go_term in go_terms_to_check:
                    if never_in_violation_found:
                        break  # Only report one conflict per constraint type
                    if check_go_term in never_in_taxon_constraints:
                        forbidden_taxa = never_in_taxon_constraints[check_go_term]
                        for forbidden_taxon in forbidden_taxa:
                            normalized_forbidden = normalize_taxon_id(forbidden_taxon)
                            # Check if given taxon or any ancestor matches the forbidden taxon
                            if normalized_forbidden in organism_ancestors:
                                inconsistent_organisms.add(organism_id)
                                inconsistencies.append((
                                    protein_id,
                                    go_term,
                                    normalized_taxon,
                                    'never_in_taxon',
                                    normalized_forbidden,
                                ))
                                never_in_violation_found = True
                                break

                # Check only_in_taxon constraints
                only_in_violation_found = False
                for check_go_term in go_terms_to_check:
                    if only_in_violation_found:
                        break  # Only report one conflict per constraint type
                    if check_go_term in in_taxon_constraints:
                        allowed_taxa = in_taxon_constraints[check_go_term]
                        # Normalize all allowed taxa
                        normalized_allowed = {normalize_taxon_id(t) for t in allowed_taxa}

                        # Check if given taxon or any ancestor is in the allowed set
                        is_allowed = False
                        conflicting_taxon = None

                        for allowed_taxon in normalized_allowed:
                            if allowed_taxon in organism_ancestors:
                                is_allowed = True
                                break

                        if not is_allowed:
                            # Report the first allowed taxon as the conflicting constraint for reference
                            conflicting_taxon = next(iter(normalized_allowed)) if normalized_allowed else "N/A"
                            inconsistent_organisms.add(organism_id)
                            inconsistencies.append((
                                protein_id,
                                go_term,
                                normalized_taxon,
                                'only_in_taxon',
                                conflicting_taxon,
                            ))
                            only_in_violation_found = True
                            break

    print(f"\nFound {len(inconsistent_organisms)} inconsistent organisms out of {len(annotations_by_organism)} total organisms")

    return inconsistencies


def save_inconsistencies(inconsistencies: List[Tuple[str, str, str, str, str]], file_path: str):
    """
    Save inconsistencies to a TSV file.
    """
    with open(file_path, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        # Write header
        writer.writerow(['protein_id', 'GO_term', 'organism_id', 'constraint_type',
                        'conflicting_taxon_id'])
        # Write data
        for inconsistency in inconsistencies:
            writer.writerow(inconsistency)


def main():
    """
    Main function to find taxon inconsistencies in test set protein annotations.
    """
    parser = argparse.ArgumentParser(
        description='Find taxon inconsistencies in test set protein annotations.'
    )
    parser.add_argument('--test-proteins', required=True,
                       help='Path to proteins_by_date_*.tsv file')
    parser.add_argument('--annotations-dir', default=DEFAULT_ANNOTATIONS_DIR,
                       help='Directory with annots_taxon_*.tsv files')
    parser.add_argument('--constraints', default=DEFAULT_CONSTRAINTS,
                       help='Path to go_taxon_constraints_extracted_obo.tsv')
    parser.add_argument('--go-hierarchy', default=DEFAULT_GO_HIERARCHY,
                       help='Path to go_hierarchy.tsv')
    parser.add_argument('--taxon-hierarchy', default=DEFAULT_TAXON_HIERARCHY,
                       help='Path to taxon_hierarchy.tsv')
    parser.add_argument('--ncbitaxon-hierarchy', default=DEFAULT_NCBITAXON_HIERARCHY,
                       help='Path to ncbitaxon_hierarchy.tsv')
    parser.add_argument('--subontology', choices=['cc', 'bp', 'mf'],
                       help='Filter by GO subontology: cc, bp, or mf')
    parser.add_argument('--go-file', default=DEFAULT_GO_FILE,
                       help='Path to go.obo file for namespace filtering')
    parser.add_argument('--output', default='test_proteins_inconsistencies.tsv',
                       help='Output file path')

    args = parser.parse_args()

    # Load GO ontology if subontology filtering is requested
    go = None
    if args.subontology:
        if not os.path.exists(args.go_file):
            print(f"[Error] GO file not found: {args.go_file}")
            print("GO file is required when using --subontology")
            sys.exit(1)

        # Add repo root to path and import Ontology from the vendored deepgo package
        if BASE_DIR not in sys.path:
            sys.path.insert(0, BASE_DIR)

        try:
            from deepgo.utils import Ontology
            print(f"Loading GO ontology from {args.go_file}...")
            go = Ontology(args.go_file)
        except ImportError as e:
            print(f"[Error] Could not import Ontology: {e}")
            print("Make sure deepgo/utils.py is available in the repo root")
            sys.exit(1)

    print("Loading test proteins...")
    test_proteins = load_test_proteins(args.test_proteins)
    print(f"Loaded {sum(len(proteins) for proteins in test_proteins.values())} test proteins across {len(test_proteins)} taxa")

    print(f"\nLoading annotations from {args.annotations_dir}...")
    annotations = load_annotations(
        args.annotations_dir,
        test_proteins,
        subontology=args.subontology,
        go=go
    )
    total_annotations = sum(
        len(go_terms)
        for proteins in annotations.values()
        for go_terms in proteins.values()
    )
    print(f"Loaded {total_annotations} annotations for {sum(len(proteins) for proteins in annotations.values())} proteins")

    print(f"\nLoading constraints from {args.constraints}...")
    in_taxon_constraints, never_in_taxon_constraints = load_constraints(args.constraints)
    print(f"Loaded {len(in_taxon_constraints)} 'only_in_taxon' constraints")
    print(f"Loaded {len(never_in_taxon_constraints)} 'never_in_taxon' constraints")

    print(f"\nLoading GO hierarchy from {args.go_hierarchy}...")
    go_hierarchy = load_go_hierarchy(args.go_hierarchy)
    print(f"Loaded {len(go_hierarchy)} GO hierarchy relationships")

    print(f"\nLoading taxon hierarchy from {args.taxon_hierarchy}...")
    taxon_is_a_hierarchy, _, _ = load_taxon_hierarchy(args.taxon_hierarchy, args.ncbitaxon_hierarchy)
    print(f"Loaded {len(taxon_is_a_hierarchy)} taxon hierarchy relationships")

    print(f"\nChecking consistency for all proteins...")
    inconsistencies = check_consistency(
        annotations,
        in_taxon_constraints,
        never_in_taxon_constraints,
        taxon_is_a_hierarchy,
        go_hierarchy
    )

    print(f"\nFound {len(inconsistencies)} inconsistencies out of {total_annotations} total GO terms")

    if inconsistencies:
        print(f"\nSaving inconsistencies to {args.output}...")
        save_inconsistencies(inconsistencies, args.output)
        print("Done!")

        # Print summary
        never_in_count = sum(1 for inc in inconsistencies if inc[3] == 'never_in_taxon')
        only_in_count = sum(1 for inc in inconsistencies if inc[3] == 'only_in_taxon')
        affected_organisms = len(set(inc[2] for inc in inconsistencies))

        print(f"\nSummary:")
        print(f"  - Total inconsistencies: {len(inconsistencies)}")
        print(f"  - never_in_taxon violations: {never_in_count}")
        print(f"  - only_in_taxon violations: {only_in_count}")
        print(f"  - Affected organisms: {affected_organisms}")
    else:
        print("\nAll annotations are consistent with taxon constraints!")


if __name__ == '__main__':
    main()
