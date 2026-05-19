#!/usr/bin/env python3
"""
GO Complex Classifier for Multiple Genomes

This script analyzes protein GO annotations for all genomes in a directory 
with a specified file extension. It counts coherent complexes (including homodimers) 
and incoherent complexes for each genome.

Usage: python complex_classifier.py <method> <file_extension> [homo_terms_file] [output_file]
"""

import sys
import os
import re
import glob
from collections import defaultdict

# Constants
MACROMOLECULAR_COMPLEX = "GO:0032991"
HOMODIMERIZATION = "GO:0042803"

def parse_homodimer_terms_file(file_path):
    """Parse the homo_terms.tsv file to identify homodimer terms.
    
    Args:
        file_path (str): Path to the homo_terms.tsv file
        
    Returns:
        set: Set of GO terms marked as homodimers ('h' or 'a')
    """
    homodimer_terms = set()
    
    try:
        with open(file_path, 'r') as f:
            # Skip header if present
            first_line = f.readline()
            if 'GO_term' in first_line and 'classification' in first_line:
                # This was a header line, continue reading
                pass
            else:
                # Not a header, process this line
                parts = first_line.strip().split()
                if len(parts) >= 2 and parts[0].startswith('GO:') and parts[1] in ['h', 'a']:
                    homodimer_terms.add(parts[0])
            
            # Process remaining lines
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0].startswith('GO:') and parts[1] in ['h', 'a']:
                    homodimer_terms.add(parts[0])
    except FileNotFoundError:
        print(f"Warning: Homodimer terms file '{file_path}' not found.")
    except Exception as e:
        print(f"Warning: Error parsing homodimer terms file: {e}")
    
    return homodimer_terms

def parse_go_ontology(file_path, include_part_of: bool = True):
    """Parse the GO ontology file to build the hierarchical structure.
    
    Args:
        file_path (str): Path to the GO ontology file (OBO format)
        
    Returns:
        dict: A dictionary mapping GO terms to their child terms
        dict: A dictionary mapping GO terms to their parent terms
        dict: A dictionary mapping GO terms to their definitions
    """
    # Initialize data structures
    term_to_children = defaultdict(set)
    term_to_parents = defaultdict(set)
    term_to_definition = {}
    
    current_term = None
    term_id = None
    
    # Parse the OBO file
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Start of a term definition
            if line == '[Term]':
                current_term = True
                term_id = None
                continue
            
            # End of term or file
            if line == '' and current_term:
                current_term = False
                continue
            
            # Skip if not in a term section
            if not current_term:
                continue
            
            # Get term ID
            if line.startswith('id:'):
                term_id = line.split('id:')[1].strip()
                continue
            
            # Get term definition
            if line.startswith('def:'):
                if term_id:
                    definition = line.split('def:')[1].strip()
                    # Extract text between quotes
                    match = re.search(r'"([^"]*)"', definition)
                    if match:
                        term_to_definition[term_id] = match.group(1)
                    else:
                        term_to_definition[term_id] = definition
                continue
            
            # Get is_a relationships
            if line.startswith('is_a:'):
                if term_id:
                    parent_id = line.split('is_a:')[1].split('!')[0].strip()
                    term_to_parents[term_id].add(parent_id)
                    term_to_children[parent_id].add(term_id)
                continue
            
            # Get part_of relationships
            if include_part_of and line.startswith('relationship: part_of'):
                if term_id:
                    parent_id = line.split('part_of')[1].split('!')[0].strip()
                    term_to_parents[term_id].add(parent_id)
                    term_to_children[parent_id].add(term_id)
    
    return term_to_children, term_to_parents, term_to_definition

def get_all_child_terms(term_id, term_to_children, recursive=True):
    """Get all child terms of a given GO term.
    
    Args:
        term_id (str): The GO term ID
        term_to_children (dict): Mapping of terms to their child terms
        recursive (bool): Whether to include all descendants (True) or only direct children (False)
        
    Returns:
        set: Set of child GO term IDs
    """
    if not recursive:
        return term_to_children.get(term_id, set())
    
    all_children = set()
    to_process = list(term_to_children.get(term_id, set()))
    
    while to_process:
        child = to_process.pop()
        if child not in all_children:
            all_children.add(child)
            to_process.extend(term_to_children.get(child, set()))
    
    return all_children

def parse_input_file(file_path):
    """Parse the input file containing protein GO term annotations.
    
    Args:
        file_path (str): Path to the input file
        
    Returns:
        dict: Mapping of protein IDs to their GO term sets
    """
    protein_go_terms = {}
    
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                protein_id = parts[0]
                go_terms = set(parts[1:])
                protein_go_terms[protein_id] = go_terms
    
    return protein_go_terms

def classify_complexes(protein_go_terms, complex_child_terms, homodimer_terms=None):
    """Classify protein complexes based on coherence.
    
    Args:
        protein_go_terms (dict): Mapping of protein IDs to GO term sets
        complex_child_terms (set): Set of child terms of GO:0032991
        homodimer_terms (set, optional): Set of GO terms known to be homodimers from external file
        
    Returns:
        dict: Mapping of complex terms to their classification (coherent/incoherent)
        dict: Mapping of complex terms to the proteins that have them
    """
    # Map of complex terms to proteins that have them
    complex_term_to_proteins = defaultdict(set)
    
    # First, identify proteins with GO:0032991 and map their complex terms
    for protein_id, terms in protein_go_terms.items():
        if MACROMOLECULAR_COMPLEX in terms:
            # Find complex child terms in this protein's annotations
            protein_complex_terms = complex_child_terms.intersection(terms)
            
            # Map each complex term to this protein
            for term in protein_complex_terms:
                complex_term_to_proteins[term].add(protein_id)
    
    # Classify each complex term
    complex_classifications = {}
    for term, proteins in complex_term_to_proteins.items():
        # Check if this term is in the list of known homodimer terms
        if homodimer_terms and term in homodimer_terms:
            # Homodimers are also considered coherent
            complex_classifications[term] = "coherent"
        # Not a homodimer, classify based on number of proteins
        elif len(proteins) > 1:
            complex_classifications[term] = "coherent"
        else:
            complex_classifications[term] = "incoherent"
    
    return complex_classifications, complex_term_to_proteins

def count_complexes(complex_classifications):
    """Count the number of coherent and incoherent complexes.
    
    Args:
        complex_classifications (dict): Mapping of complex terms to their classification
        
    Returns:
        tuple: (coherent_count, incoherent_count)
    """
    coherent_count = sum(1 for c in complex_classifications.values() if c == "coherent")
    incoherent_count = sum(1 for c in complex_classifications.values() if c == "incoherent")
    
    return coherent_count, incoherent_count

def extract_genome_name(file_path, file_extension):
    """Extract the genome name from a file path.
    
    Args:
        file_path (str): Full path to the genome file
        file_extension (str): File extension to remove
        
    Returns:
        str: Genome name
    """
    # Get the base filename without directory
    base_name = os.path.basename(file_path)
    genome_name = base_name
    
    # Remove "_IPscan_GO_ancestors" if present
    if genome_name.endswith("_deepgometa_th_GO_ancestors.tsv"):
        genome_name = genome_name[:-len("_deepgometa_th_GO_ancestors.tsv")]
    
    return genome_name

def process_genome_file(file_path, method, file_extension, complex_child_terms, homodimer_terms=None):
    """Process a single genome file.
    
    Args:
        file_path (str): Path to the genome file
        method (str): Method name
        file_extension (str): File extension
        complex_child_terms (set): Set of child terms of GO:0032991
        homodimer_terms (set, optional): Set of GO terms known to be homodimers
        
    Returns:
        tuple: (genome, method, incoherent_count, coherent_count)
    """
    try:
        # Extract genome name from file path
        genome = extract_genome_name(file_path, file_extension)
        
        # Parse protein annotations
        protein_go_terms = parse_input_file(file_path)
        
        # Classify complexes
        complex_classifications, _ = classify_complexes(
            protein_go_terms, complex_child_terms, homodimer_terms
        )
        
        # Count complex types
        coherent_count, incoherent_count = count_complexes(complex_classifications)
        
        return genome, method, incoherent_count, coherent_count
    except Exception as e:
        print(f"Error processing file {file_path}: {e}", file=sys.stderr)
        genome = extract_genome_name(file_path, file_extension)
        return genome, method, 0, 0

def main():
    """Main entry point for the script."""
    if len(sys.argv) < 3:
        print(f"Usage: {os.path.basename(sys.argv[0])} <method> <file_extension> [homo_terms_file] [output_file]")
        sys.exit(1)
    
    method = sys.argv[1]
    file_extension = sys.argv[2]
    homo_terms_file = None
    output_file = f"{method}_complex_results.tsv"  # Default output filename
    
    # Parse command line arguments
    if len(sys.argv) > 3:
        homo_terms_file = sys.argv[3]
    
    if len(sys.argv) > 4:
        output_file = sys.argv[4]
    
    # Ensure file extension starts with a dot if needed
    if not file_extension.startswith('.') and not file_extension.startswith('*'):
        file_extension = file_extension
    
    # Directory where genome files are located
    #genome_dir = f"../preds_{method}/new_GO/"
    genome_dir = f"../model_organisms/IEA/"
    
    # Check if directory exists
    if not os.path.exists(genome_dir):
        print(f"Error: Directory '{genome_dir}' not found.")
        sys.exit(1)
    
    # Default ontology file path
    ontology_file = '../data/go-basic.obo'
    
    # Check if ontology file exists
    if not os.path.exists(ontology_file):
        print(f"Error: Ontology file '{ontology_file}' not found.")
        sys.exit(1)
    
    # Parse homo_terms file if provided
    homodimer_terms = None
    if homo_terms_file:
        if not os.path.exists(homo_terms_file):
            print(f"Warning: Homodimer terms file '{homo_terms_file}' not found.")
        else:
            print(f"Parsing homodimer terms file: {homo_terms_file}")
            homodimer_terms = parse_homodimer_terms_file(homo_terms_file)
            print(f"Found {len(homodimer_terms)} known homodimer terms")
    
    # Parse GO ontology to get term relationships
    print(f"Parsing GO ontology file: {ontology_file}")
    term_to_children, _, _ = parse_go_ontology(ontology_file, include_part_of=False)
    
    # Get all child terms of GO:0032991 (macromolecular complex)
    print(f"Identifying child terms of {MACROMOLECULAR_COMPLEX} (macromolecular complex)")
    complex_child_terms = get_all_child_terms(MACROMOLECULAR_COMPLEX, term_to_children)
    complex_child_terms.add(MACROMOLECULAR_COMPLEX)  # Include the parent term itself
    print(f"Found {len(complex_child_terms)} child terms")
    
    # Find all genome files with the specified extension
    genome_pattern = os.path.join(genome_dir, f"*{file_extension}")
    genome_files = glob.glob(genome_pattern)
    
    if not genome_files:
        print(f"Error: No files found matching pattern '{genome_pattern}'.")
        sys.exit(1)
    
    print(f"Found {len(genome_files)} genome files to process")
    
    # Process each genome file and collect results
    results = []
    for file_path in genome_files:
        print(f"Processing {os.path.basename(file_path)}...")
        result = process_genome_file(file_path, method, file_extension, complex_child_terms, homodimer_terms)
        results.append(result)
    
    # Write results to output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("genome\tmethod\tincoherent_complexes\tcoherent_complexes\n")
        
        # Write each result row
        for genome, method, incoherent, coherent in sorted(results):
            f.write(f"{genome}\t{method}\t{incoherent}\t{coherent}\n")
    
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    main()