"""
Parse the ncbitaxon.obo file and output the hierarchy to a tsv file.
"""

import argparse
import os
import re
from collections import defaultdict

def parse_ncbitaxon(ncbitaxon_obo_file):
    """
    Parse the ncbitaxon.obo file and output the hierarchy to a tsv file.
    """
    hierarchy = defaultdict(set)
    with open(ncbitaxon_obo_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('id:'):
                current_term = line.split('id:')[1].strip()
                current_term = current_term.replace(":", "_")
            elif line.startswith('is_a:'):
                parent_term = line.split('is_a:')[1].split('!')[0].strip()
                parent_term = parent_term.replace(":", "_")
                hierarchy[current_term].add(parent_term)
    return hierarchy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse ncbitaxon.obo file and output the hierarchy to a tsv file')
    parser.add_argument('--ncbitaxon-obo-file', type=str, default=None, help='Path to the ncbitaxon.obo file')
    parser.add_argument('--output-file', type=str, default=None, help='Path to the output file')
    args = parser.parse_args()

    if args.ncbitaxon_obo_file is None:
        args.ncbitaxon_obo_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'ncbitaxon.obo')
    if args.output_file is None:
        args.output_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'ncbitaxon_hierarchy.tsv')

    hierarchy = parse_ncbitaxon(args.ncbitaxon_obo_file)
    print(f"[INFO] Parsed {len(hierarchy)} ncbitaxon terms and {sum(len(parents) for parents in hierarchy.values())} parent-child relationships")
    with open(args.output_file, 'w') as f:
        for term in hierarchy:
            for parent in hierarchy[term]:
                f.write(f"{term}\t{parent}\n")
    print(f"[INFO] Saved to {args.output_file}")

