"""
Given a GO OBO file, parse it and output the hierarchy in a format that can be used by the problem formulation.
The hierarchy file is in the format:
GO:0005737	GO:0005623
GO:0005737	GO:0005624

where the first column is the child term and the second column is the parent term for is_a relation.

TODO: expand to part_of, has_part, and other relations.
"""

import argparse
import os
import re
from collections import defaultdict

NAMESPACES = {
    'mf': 'molecular_function',
    'bp': 'biological_process',
    'cc': 'cellular_component',
}

def parse_go_hierarchy(go_obo_file, subontology):
    """
    Parse the GO OBO file and return the hierarchy in a format that can be used by the problem formulation.
    """
    if subontology is not None and subontology not in NAMESPACES:
        raise ValueError(f"Invalid subontology: {subontology}. Valid subontologies are: {NAMESPACES.keys()}")
    if subontology is None:
        namespace = None
    else:
        namespace = NAMESPACES[subontology]
    hierarchy = defaultdict(set)
    with open(go_obo_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('id:'):
                current_term = line.split('id:')[1].strip()
                to_be_added = True
            elif line.startswith('namespace:'):
                line_namespace = line.split('namespace:')[1].strip()
                if namespace is not None and line_namespace != namespace:
                    to_be_added = False
            elif line.startswith('is_a:'):
                parent_term = line.split('is_a:')[1].split('!')[0].strip()
                if to_be_added:
                    hierarchy[current_term].add(parent_term)
    return hierarchy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse GO OBO file and output the hierarchy')
    parser.add_argument('--go-obo-file', type=str, default=None, help='Path to the GO OBO file')
    parser.add_argument('--output-file', type=str, default=None, help='Path to the output file')
    parser.add_argument('--subontology', type=str, default=None, help='Subontology to parse (mf, bp, cc)')
    args = parser.parse_args()

    if args.go_obo_file is None:
        args.go_obo_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'go-basic.obo')
    if args.output_file is None:
        if args.subontology is None:
            args.output_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'go_hierarchy.tsv')
        else:
            args.output_file = os.path.join(os.path.dirname(__file__), '..', 'data', f'go_hierarchy_{args.subontology}.tsv')

    hierarchy = parse_go_hierarchy(args.go_obo_file, args.subontology)
    print(f"[INFO] Parsed {len(hierarchy)} GO terms and {sum(len(parents) for parents in hierarchy.values())} parent-child relationships")
    with open(args.output_file, 'w') as f:
        for term in hierarchy:
            for parent in hierarchy[term]:
                f.write(f"{term}\t{parent}\n")
    print(f"[INFO] Saved to {args.output_file}")