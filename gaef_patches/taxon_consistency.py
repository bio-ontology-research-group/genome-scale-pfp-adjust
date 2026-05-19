#!/usr/bin/env python3
"""
Taxon Consistency Checker - Python implementation of taxon_consistency.groovy

This module checks whether GO annotation taxon constraints are satisfiable
for a given genome. It replaces the Groovy/OWLAPI/ELK reasoner approach
with custom Python logic for better performance.

Usage:
    python taxon_consistency.py <input_file> <output_file>
"""

import sys
import os
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set, Tuple, List, NamedTuple, Optional
from dataclasses import dataclass, field


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class GenomeAnnotation:
    """Represents a single protein annotation with taxon constraints."""
    genome_name: str
    protein_name: str
    go_id: str
    never_in_taxon: Set[str] = field(default_factory=set)
    only_in_taxon: Set[str] = field(default_factory=set)


@dataclass
class ConsistencyResult:
    """Result of a satisfiability check."""
    is_satisfiable: bool
    explanation: str = ""


# ============================================================================
# OWL Parsing (Phase 1)
# ============================================================================

def parse_owl_files(taxon_file: str, go_taxon_file: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]]]:
    """
    Parse OWL files to extract taxonomy hierarchy and disjointness axioms.
    
    This is a lightweight XML parser optimized for the specific OWL patterns
    used in these files, avoiding the overhead of full rdflib parsing.
    
    Returns:
        Tuple of (subclass_of, disjoint_with, union_members) dictionaries
    """
    import xml.etree.ElementTree as ET
    
    # Namespaces used in the OWL files
    ns = {
        'owl': 'http://www.w3.org/2002/07/owl#',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
    }
    
    OBO_PREFIX = "http://purl.obolibrary.org/obo/"
    
    # Data structures
    subclass_of: Dict[str, Set[str]] = defaultdict(set)  # child -> parents
    disjoint_with: Dict[str, Set[str]] = defaultdict(set)  # taxon -> disjoint taxa
    union_members: Dict[str, Set[str]] = defaultdict(set)  # union class -> member taxa
    
    def extract_taxon_id(uri: str) -> Optional[str]:
        """Extract taxon ID from URI."""
        if uri.startswith(OBO_PREFIX):
            return uri[len(OBO_PREFIX):]
        return None
    
    def parse_file(filepath: str):
        """Parse a single OWL file."""
        tree = ET.parse(filepath)
        root = tree.getroot()
        
        # Process all elements
        for elem in root.iter():
            tag = elem.tag.split('}')[-1] if '}' in elem.tag else elem.tag
            
            # Handle rdf:Description elements (disjointness in ncbitaxon file)
            if tag == 'Description':
                about = elem.get(f"{{{ns['rdf']}}}about")
                if about:
                    subject_id = extract_taxon_id(about)
                    if subject_id:
                        for child in elem:
                            child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                            if child_tag == 'disjointWith':
                                resource = child.get(f"{{{ns['rdf']}}}resource")
                                if resource:
                                    object_id = extract_taxon_id(resource)
                                    if object_id:
                                        disjoint_with[subject_id].add(object_id)
                                        disjoint_with[object_id].add(subject_id)  # symmetric
            
            # Handle owl:Class elements
            elif tag == 'Class':
                about = elem.get(f"{{{ns['rdf']}}}about")
                if about:
                    subject_id = extract_taxon_id(about)
                    if subject_id:
                        for child in elem:
                            child_tag = child.tag.split('}')[-1] if '}' in child.tag else child.tag
                            
                            # rdfs:subClassOf
                            if child_tag == 'subClassOf':
                                resource = child.get(f"{{{ns['rdf']}}}resource")
                                if resource:
                                    parent_id = extract_taxon_id(resource)
                                    if parent_id:
                                        subclass_of[subject_id].add(parent_id)
                            
                            # owl:disjointWith
                            elif child_tag == 'disjointWith':
                                resource = child.get(f"{{{ns['rdf']}}}resource")
                                if resource:
                                    object_id = extract_taxon_id(resource)
                                    if object_id:
                                        disjoint_with[subject_id].add(object_id)
                                        disjoint_with[object_id].add(subject_id)  # symmetric
                            
                            # owl:equivalentClass (for union definitions)
                            elif child_tag == 'equivalentClass':
                                # Look for owl:unionOf inside
                                for desc in child.iter():
                                    desc_tag = desc.tag.split('}')[-1] if '}' in desc.tag else desc.tag
                                    if desc_tag == 'unionOf':
                                        for member in desc:
                                            member_about = member.get(f"{{{ns['rdf']}}}about")
                                            if member_about:
                                                member_id = extract_taxon_id(member_about)
                                                if member_id:
                                                    union_members[subject_id].add(member_id)
    
    # Parse both files
    if os.path.exists(taxon_file):
        parse_file(taxon_file)
    if os.path.exists(go_taxon_file):
        parse_file(go_taxon_file)
    
    return dict(subclass_of), dict(disjoint_with), dict(union_members)


# ============================================================================
# Taxon Hierarchy (Phase 2)
# ============================================================================

class TaxonHierarchy:
    """
    Manages taxon hierarchy and disjointness relationships.
    Provides efficient lookups for satisfiability checking.
    """
    
    def __init__(self, subclass_of: Dict[str, Set[str]], 
                 disjoint_with: Dict[str, Set[str]], 
                 union_members: Dict[str, Set[str]]):
        self.subclass_of = subclass_of
        self.disjoint_with = disjoint_with
        self.union_members = union_members
        
        # Pre-compute transitive superclasses for each taxon
        self._ancestors_cache: Dict[str, Set[str]] = {}
        
        # Build reverse mapping (parent -> children)
        self.subclass_from: Dict[str, Set[str]] = defaultdict(set)
        for child, parents in subclass_of.items():
            for parent in parents:
                self.subclass_from[parent].add(child)
    
    def get_ancestors(self, taxon_id: str) -> Set[str]:
        """Get all ancestors (superclasses) of a taxon, including itself."""
        if taxon_id in self._ancestors_cache:
            return self._ancestors_cache[taxon_id]
        
        ancestors = {taxon_id}
        to_visit = list(self.subclass_of.get(taxon_id, set()))
        
        while to_visit:
            current = to_visit.pop()
            if current not in ancestors:
                ancestors.add(current)
                to_visit.extend(self.subclass_of.get(current, set()))
        
        self._ancestors_cache[taxon_id] = ancestors
        return ancestors
    
    def get_descendants(self, taxon_id: str) -> Set[str]:
        """Get all descendants (subclasses) of a taxon, including itself."""
        descendants = {taxon_id}
        to_visit = list(self.subclass_from.get(taxon_id, set()))
        
        while to_visit:
            current = to_visit.pop()
            if current not in descendants:
                descendants.add(current)
                to_visit.extend(self.subclass_from.get(current, set()))
        
        return descendants
    
    def expand_union(self, taxon_id: str) -> Set[str]:
        """Expand a union class to its constituent taxa."""
        if taxon_id in self.union_members:
            result = set()
            for member in self.union_members[taxon_id]:
                result.update(self.expand_union(member))
            return result
        return {taxon_id}
    
    def is_subclass_of(self, taxon_a: str, taxon_b: str) -> bool:
        """Check if taxon_a is a subclass of taxon_b."""
        return taxon_b in self.get_ancestors(taxon_a)
    
    def are_disjoint(self, taxon_a: str, taxon_b: str) -> Tuple[bool, Optional[Tuple[str, str]]]:
        """
        Check if two taxa are disjoint.
        Returns (is_disjoint, (taxon1, taxon2) if disjoint axiom found else None)
        """
        # Direct disjointness
        if taxon_b in self.disjoint_with.get(taxon_a, set()):
            return True, (taxon_a, taxon_b)
        
        # Check if any ancestor is disjoint
        ancestors_a = self.get_ancestors(taxon_a)
        ancestors_b = self.get_ancestors(taxon_b)
        
        for anc_a in ancestors_a:
            for disjoint in self.disjoint_with.get(anc_a, set()):
                if disjoint in ancestors_b:
                    return True, (anc_a, disjoint)
        
        for anc_b in ancestors_b:
            for disjoint in self.disjoint_with.get(anc_b, set()):
                if disjoint in ancestors_a:
                    return True, (anc_b, disjoint)
        
        return False, None
    
    def normalize_taxon_id(self, taxon_id: str) -> str:
        """Normalize taxon ID by removing common prefixes."""
        return taxon_id.replace("NCBITaxon_", "").replace("NCBITaxon_Union_", "Union_")


# ============================================================================
# Satisfiability Checking (Phase 3)
# ============================================================================

def check_genome_satisfiability(
    annotations: List[GenomeAnnotation],
    hierarchy: TaxonHierarchy
) -> ConsistencyResult:
    """
    Check if the taxon constraints for a genome are satisfiable.
    
    The Groovy script creates an OWL intersection of:
    - All "only_in_taxon" classes (organism must be in ALL of these)
    - All negated "never_in_taxon" classes (using _neg classes - organism must NOT be in these)
    
    The intersection is unsatisfiable if:
    1. Two "only_in" taxa are disjoint (can't be in both)
    2. An "only_in" taxon equals or is a subclass of a "never_in" taxon
    3. An "only_in" taxon's _neg class would need to be in the intersection
       but conflicts with the "only" requirement
       
    The key insight: if we have "only_in X" and "never_in X", 
    the intersection would require being in X AND in X_neg, 
    which is impossible since X and X_neg are defined as disjoint.
    """
    
    # Collect all constraints
    only_set: Set[str] = set()
    never_set: Set[str] = set()
    
    # Also track which annotations contribute which constraints
    only_annotations: Dict[str, List[GenomeAnnotation]] = defaultdict(list)
    never_annotations: Dict[str, List[GenomeAnnotation]] = defaultdict(list)
    
    for ann in annotations:
        for taxon in ann.only_in_taxon:
            normalized = normalize_taxon_id(taxon)
            only_set.add(normalized)
            only_annotations[normalized].append(ann)
        
        for taxon in ann.never_in_taxon:
            normalized = normalize_taxon_id(taxon)
            never_set.add(normalized)
            never_annotations[normalized].append(ann)
    
    # If no constraints, trivially satisfiable
    if not only_set and not never_set:
        return ConsistencyResult(is_satisfiable=True)
    
    # Check 1: Direct conflict - same taxon in both "only" and "never"
    # If we require "only_in X" and "never_in X", we need X ∩ X_neg which is empty
    common_taxa = only_set & never_set
    if common_taxa:
        taxon = next(iter(common_taxa))
        explanation = generate_explanation(
            conflict_type="only_equals_never",
            taxon_a=taxon, taxon_b=taxon,
            only_annotations=only_annotations,
            never_annotations=never_annotations,
            disjoint_pair=(f"NCBITaxon_{taxon}", f"NCBITaxon_{taxon}")
        )
        return ConsistencyResult(is_satisfiable=False, explanation=explanation)
    
    # Check 2: Are all "only_in" constraints mutually compatible?
    # For the intersection to be non-empty, the taxa must not be disjoint
    only_list = list(only_set)
    for i, taxon_a in enumerate(only_list):
        full_a = make_full_taxon_id(taxon_a)
        for taxon_b in only_list[i+1:]:
            full_b = make_full_taxon_id(taxon_b)
            is_disjoint, disjoint_pair = hierarchy.are_disjoint(full_a, full_b)
            if is_disjoint:
                explanation = generate_explanation(
                    conflict_type="only_only_disjoint",
                    taxon_a=taxon_a, taxon_b=taxon_b,
                    only_annotations=only_annotations,
                    never_annotations=never_annotations,
                    disjoint_pair=disjoint_pair
                )
                return ConsistencyResult(is_satisfiable=False, explanation=explanation)
    
    # Check 3: Do "never_in" constraints conflict with "only_in" constraints?
    # If an "only" taxon is a subclass of (or equal to) a "never" taxon's ancestor,
    # then being in "only" means being in "never" which is forbidden
    for only_taxon in only_set:
        full_only = make_full_taxon_id(only_taxon)
        only_ancestors = hierarchy.get_ancestors(full_only)
        
        for never_taxon in never_set:
            full_never = make_full_taxon_id(never_taxon)
            never_ancestors = hierarchy.get_ancestors(full_never)
            
            # If only_taxon is a subclass of never_taxon, conflict
            # (being in only means being in never, but never is forbidden)
            if full_never in only_ancestors:
                explanation = generate_explanation(
                    conflict_type="only_subclass_of_never",
                    taxon_a=only_taxon, taxon_b=never_taxon,
                    only_annotations=only_annotations,
                    never_annotations=never_annotations,
                    disjoint_pair=(full_only, full_never)
                )
                return ConsistencyResult(is_satisfiable=False, explanation=explanation)
            
            # If never_taxon is a subclass of only_taxon, we need to check
            # if all organisms in only_taxon are also in never_taxon
            # This happens when they're the same or when only is more general
            if full_only in never_ancestors:
                # Only if only_taxon == never_taxon (already checked above)
                # or if there's no room for organisms that are in only but not in never
                # For now, this is handled by the disjointness check below
                pass
    
    # Check 4: Check if any "only" taxon is disjoint with any "never" taxon's negation
    # The intersection needs: only_taxon AND never_taxon_neg
    # If only_taxon is disjoint with never_taxon_neg, the intersection is empty
    for only_taxon in only_set:
        full_only = make_full_taxon_id(only_taxon)
        
        for never_taxon in never_set:
            neg_class = f"NCBITaxon_{never_taxon}_neg"
            
            # Check if only taxon is disjoint with the neg class
            is_disjoint, disjoint_pair = hierarchy.are_disjoint(full_only, neg_class)
            if is_disjoint:
                explanation = generate_explanation(
                    conflict_type="only_disjoint_with_never_neg",
                    taxon_a=only_taxon, taxon_b=never_taxon,
                    only_annotations=only_annotations,
                    never_annotations=never_annotations,
                    disjoint_pair=disjoint_pair
                )
                return ConsistencyResult(is_satisfiable=False, explanation=explanation)
    
    return ConsistencyResult(is_satisfiable=True)


def normalize_taxon_id(taxon_id: str) -> str:
    """Normalize taxon ID by removing NCBITaxon_ prefix."""
    return taxon_id.replace("NCBITaxon_", "")


def make_full_taxon_id(normalized_id: str) -> str:
    """Convert normalized ID back to full IRI fragment."""
    if normalized_id.startswith("Union_"):
        return f"NCBITaxon_{normalized_id}"
    return f"NCBITaxon_{normalized_id}"


# ============================================================================
# Explanation Generation (Phase 4)
# ============================================================================

def generate_explanation(
    conflict_type: str,
    taxon_a: str,
    taxon_b: str,
    only_annotations: Dict[str, List[GenomeAnnotation]],
    never_annotations: Dict[str, List[GenomeAnnotation]],
    disjoint_pair: Optional[Tuple[str, str]] = None
) -> str:
    """
    Generate a human-readable explanation for an unsatisfiable result.
    Matches the format of the Groovy script output.
    """
    parts = []
    
    # Collect annotations that contribute to the conflict
    if taxon_a in only_annotations:
        annots = only_annotations[taxon_a]
        protein_list = ", ".join(f"Protein {a.protein_name} ({a.go_id})" for a in annots[:50])
        if len(annots) > 50:
            protein_list += f", ... and {len(annots) - 50} more"
        parts.append(f"{protein_list} requires only in taxon T{taxon_a}")
    
    if taxon_b in never_annotations:
        annots = never_annotations[taxon_b]
        protein_list = ", ".join(f"Protein {a.protein_name} ({a.go_id})" for a in annots[:50])
        if len(annots) > 50:
            protein_list += f", ... and {len(annots) - 50} more"
        parts.append(f"{protein_list} requires never in taxon T{taxon_b}")
    
    if taxon_a in never_annotations:
        annots = never_annotations[taxon_a]
        protein_list = ", ".join(f"Protein {a.protein_name} ({a.go_id})" for a in annots[:50])
        if len(annots) > 50:
            protein_list += f", ... and {len(annots) - 50} more"
        parts.append(f"{protein_list} requires never in taxon T{taxon_a}")
    
    if taxon_b in only_annotations:
        annots = only_annotations[taxon_b]
        protein_list = ", ".join(f"Protein {a.protein_name} ({a.go_id})" for a in annots[:50])
        if len(annots) > 50:
            protein_list += f", ... and {len(annots) - 50} more"
        parts.append(f"{protein_list} requires only in taxon T{taxon_b}")
    
    explanation = "; ".join(parts)
    
    # Add disjointness relationship if available
    if disjoint_pair:
        id1 = disjoint_pair[0].replace("NCBITaxon_", "").replace("_neg", "")
        id2 = disjoint_pair[1].replace("NCBITaxon_", "").replace("_neg", "")
        explanation += f", and T{id1} and T{id2} are disjoint"
    
    return explanation


# ============================================================================
# Input Parsing
# ============================================================================

def parse_input_file(input_file: str) -> Dict[str, List[GenomeAnnotation]]:
    """
    Parse the input TSV file and group annotations by genome.
    
    Input format:
        protein_name	GO_ID	never_in_taxon	only_in_taxon
    
    The genome name is derived from the input filename.
    """
    genome_annotations: Dict[str, List[GenomeAnnotation]] = defaultdict(list)
    
    # Derive genome name from filename
    genome_name = Path(input_file).stem.replace("_consistency", "")
    
    with open(input_file, 'r', encoding='utf-8') as f:
        # Skip header
        header = f.readline()
        
        for line in f:
            # Use split with -1 limit equivalent to preserve trailing empty fields
            # Or pad the parts list if needed
            parts = line.rstrip('\n\r').split('\t')
            
            # Need at least protein_name and go_id
            if len(parts) < 2:
                continue
            
            protein_name = parts[0].strip()
            go_id = parts[1].strip()
            
            # Handle cases where columns might be missing
            never_col = parts[2] if len(parts) > 2 else ""
            only_col = parts[3] if len(parts) > 3 else ""
            
            never_taxa = set(t.strip() for t in never_col.split(',') if t.strip())
            only_taxa = set(t.strip() for t in only_col.split(',') if t.strip())
            
            # Only include annotations that have at least one constraint
            if never_taxa or only_taxa:
                annotation = GenomeAnnotation(
                    genome_name=genome_name,
                    protein_name=protein_name,
                    go_id=go_id,
                    never_in_taxon=never_taxa,
                    only_in_taxon=only_taxa
                )
                
                genome_annotations[genome_name].append(annotation)
    
    return genome_annotations


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main entry point for command-line usage."""
    if len(sys.argv) != 3:
        print("Usage: python taxon_consistency.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Vendored constraints live at <repo>/data/constraints/.
    constraints_dir = Path(__file__).resolve().parent.parent / "data" / "constraints"
    taxon_file = constraints_dir / "ncbitaxon_with_disjointness.owl"
    go_taxon_file = constraints_dir / "go-taxon-groupings.owl"
    
    # Verify input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)
    
    # Verify ontology files exist
    if not taxon_file.exists() or not go_taxon_file.exists():
        print("Error: Ontology file(s) not found")
        print(f"  Expected: {taxon_file}")
        print(f"  Expected: {go_taxon_file}")
        sys.exit(1)
    
    print("Loading ontologies...")
    subclass_of, disjoint_with, union_members = parse_owl_files(
        str(taxon_file), str(go_taxon_file)
    )
    print(f"  Loaded {len(subclass_of)} subclass relationships")
    print(f"  Loaded {len(disjoint_with)} disjointness axioms")
    print(f"  Loaded {len(union_members)} union class definitions")
    
    # Build hierarchy
    hierarchy = TaxonHierarchy(subclass_of, disjoint_with, union_members)
    print("Ontologies loaded and processed successfully.")
    
    # Parse input file
    genome_annotations = parse_input_file(input_file)
    print(f"Found {len(genome_annotations)} genome(s).")
    
    # Write output header
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Genome\tIsSatisfiable\tExplanation\n")
    
    # Process each genome
    for genome_name, annotations in genome_annotations.items():
        print(f"\nProcessing genome: {genome_name}")
        
        result = check_genome_satisfiability(annotations, hierarchy)
        
        # Write result
        with open(output_file, 'a', encoding='utf-8') as f:
            explanation = result.explanation.replace('\n', ' ').replace('\t', ' ')
            f.write(f"{genome_name}\t{result.is_satisfiable}\t{explanation}\n")
    
    print("\nFinished processing genomes.")


# ============================================================================
# API for programmatic use
# ============================================================================

def check_taxon_consistency(
    input_file: str,
    output_file: str,
    constraints_dir: Optional[str] = None
) -> Dict[str, ConsistencyResult]:
    """
    Check taxon consistency for annotations in input_file.
    
    This is the main API function for use from other Python code.
    
    Args:
        input_file: Path to TSV file with protein annotations
        output_file: Path for output TSV file
        constraints_dir: Directory containing OWL constraint files
                        (defaults to script's constraints/ directory)
    
    Returns:
        Dictionary mapping genome names to ConsistencyResult objects
    """
    # Determine constraint file locations (default = vendored <repo>/data/constraints/).
    if constraints_dir is None:
        constraints_dir = Path(__file__).resolve().parent.parent / "data" / "constraints"
    else:
        constraints_dir = Path(constraints_dir)
    
    taxon_file = constraints_dir / "ncbitaxon_with_disjointness.owl"
    go_taxon_file = constraints_dir / "go-taxon-groupings.owl"
    
    # Parse ontologies
    subclass_of, disjoint_with, union_members = parse_owl_files(
        str(taxon_file), str(go_taxon_file)
    )
    
    # Build hierarchy
    hierarchy = TaxonHierarchy(subclass_of, disjoint_with, union_members)
    
    # Parse input
    genome_annotations = parse_input_file(input_file)
    
    # Write output header
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Genome\tIsSatisfiable\tExplanation\n")
    
    # Process each genome
    results = {}
    for genome_name, annotations in genome_annotations.items():
        result = check_genome_satisfiability(annotations, hierarchy)
        results[genome_name] = result
        
        # Write result
        with open(output_file, 'a', encoding='utf-8') as f:
            explanation = result.explanation.replace('\n', ' ').replace('\t', ' ')
            f.write(f"{genome_name}\t{result.is_satisfiable}\t{explanation}\n")
    
    return results


if __name__ == '__main__':
    main()

