"""
Given a go-taxon-groupings.obo file, parse the taxon hierarchy and output the hierarchy in a format that can be used by the problem formulation.

TODO:
- add support for ncbitaxon.obo file
    - using only relevant taxons in relation to go_taxon_constraints_updated.tsv
"""

import argparse
import os
from collections import defaultdict

def parse_taxon_hierarchy(go_taxon_groupings_file):
    """
    Parse the taxon hierarchy from the go-taxon-groupings.obo file.
    """
    is_a_hierarchy = defaultdict(set)
    disjoint_from_hierarchy = defaultdict(set)
    union_of_hierarchy = defaultdict(set)
    with open(go_taxon_groupings_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('id:'):
                current_term = line.split('id:')[1].strip()
                current_term = current_term.replace(":", "_")
            elif line.startswith('is_a:'):
                parent_term = line.split('is_a:')[1].split('!')[0].strip()
                parent_term = parent_term.replace(":", "_")
                is_a_hierarchy[current_term].add(parent_term)
            elif line.startswith('disjoint_from:'):
                disjoint_from_term = line.split('disjoint_from:')[1].split('!')[0].strip()
                disjoint_from_term = disjoint_from_term.replace(":", "_")
                disjoint_from_hierarchy[current_term].add(disjoint_from_term)
            elif line.startswith('union_of:'):
                taxon_term = line.split('union_of:')[1].split('!')[0].strip()
                taxon_term = taxon_term.replace(":", "_")
                union_of_hierarchy[current_term].add(taxon_term)
    return is_a_hierarchy, disjoint_from_hierarchy, union_of_hierarchy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse taxon hierarchy')
    parser.add_argument('--go-taxon-groupings-file', type=str, default="data/go-taxon-groupings.obo", help='Path to the go-taxon-groupings.obo file')
    parser.add_argument('--output-file', type=str, default="data/taxon_hierarchy.tsv", help='Path to the output file')
    args = parser.parse_args()

    is_a_hierarchy, disjoint_from_hierarchy, union_of_hierarchy = parse_taxon_hierarchy(args.go_taxon_groupings_file)
    print(f"[INFO] Parsed {len(is_a_hierarchy)} is_a relationships and {len(disjoint_from_hierarchy)} disjoint_from relationships and {len(union_of_hierarchy)} union_of relationships")
    with open(args.output_file, 'w') as f:
        f.write("Term\tRelationship\tParent/Disjoint_From_Term\n")
        for term in is_a_hierarchy:
            for parent in is_a_hierarchy[term]:
                f.write(f"{term}\tis_a\t{parent}\n")
        for term in disjoint_from_hierarchy:
            for disjoint_from in disjoint_from_hierarchy[term]:
                f.write(f"{term}\tdisjoint_from\t{disjoint_from}\n")
        for term in union_of_hierarchy:
            for union_of in union_of_hierarchy[term]:
                f.write(f"{term}\tunion_of\t{union_of}\n")
    print(f"[INFO] Saved to {args.output_file}")