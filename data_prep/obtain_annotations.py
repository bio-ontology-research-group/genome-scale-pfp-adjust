"""
Given the proteomes ids, obtain propagated experimental annotations from .dat files.
"""

import gzip
import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from deepgo.utils import Ontology


EXP_CODES = set([
    'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC',
    'HTP', 'HDA', 'HMP', 'HGI', 'HEP'])


def is_exp_code(code):
    return code in EXP_CODES

def load_annotations(proteome_dat_gz_file):
    """
    Parses proteome .dat.gz file and loads list of proteins and their annotations.
    Args:
       proteome_dat_gz_file (string): A path to the data file
    Returns:
       Tuple of 2 lists (proteins, exp_annotations)
    """
    proteins = list()
    exp_annotations = list()

    with gzip.open(proteome_dat_gz_file, 'rt') as f:
        prot_id = ''
        annots = list()
        for line in f:
            items = line.strip().split('   ')
            if items[0] == 'ID' and len(items) > 1:
                if prot_id != '':
                    proteins.append(prot_id)
                    exp_annotations.append(annots)
                prot_id = items[1]
                annots = list()

            elif items[0] == 'DR' and len(items) > 1:
                items = items[1].split('; ')
                if items[0] == 'GO':
                    go_id = items[1]
                    code = items[3].split(':')[0]
                    annots.append(go_id + '|' + code)

        proteins.append(prot_id)
        exp_annotations.append(annots)

    # filter experimental annotations
    filtered_exp_annotations = list()
    for annots in exp_annotations:
        filtered_annots = list()
        for annot in annots:
            go_id, code = annot.split('|')
            if is_exp_code(code):
                filtered_annots.append(go_id)
        filtered_exp_annotations.append(filtered_annots)

    return proteins, filtered_exp_annotations

def propagate_annotations(exp_annotations, go_ontology):
    """
    Propagates experimental annotations to ancestors.
    Args:
       exp_annotations (list): List of experimental annotations
       go_ontology (Ontology): GO ontology
    Returns:
       List of propagated annotations
    """
    prop_annotations = list()
    for protein_exp_annots in exp_annotations:
        protein_prop_annots = set()
        for annot in protein_exp_annots:
            protein_prop_annots |= go_ontology.get_ancestors(annot)
        prop_annotations.append(list(protein_prop_annots))

    return prop_annotations

def save_annotations(proteins, prop_annotations, output_file):
    """
    Saves list of proteins and their experimental annotations to a TSV file.
    Args:
       proteins (list): List of proteins
       prop_annotations (list): List of propagated annotations
       output_file (string): A path to the output file
    """
    with open(output_file, 'w') as f:
        for protein, annotation in zip(proteins, prop_annotations):
            f.write(f"{protein}\t" + "\t".join(annotation) + "\n")


def main(proteome_file, output_file, go_file):
    print(f"Loading GO ontology from {go_file}")
    go_ontology = Ontology(go_file, with_rels=True)

    proteins, exp_annotations = load_annotations(proteome_file)
    prop_annotations = propagate_annotations(exp_annotations, go_ontology)

    save_annotations(proteins, prop_annotations, output_file)
    print(f"Saved {len(proteins)} proteins to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--proteome_file', type=str, required=True)
    parser.add_argument('--output_file', type=str, default=None)
    parser.add_argument('--go_file', type=str, default=None)
    args = parser.parse_args()

    if args.output_file is None:
        args.output_file = os.path.splitext(args.proteome_file)[0] + '.tsv'
    if args.go_file is None:
        args.go_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'go.obo')

    main(args.proteome_file, args.output_file, args.go_file)