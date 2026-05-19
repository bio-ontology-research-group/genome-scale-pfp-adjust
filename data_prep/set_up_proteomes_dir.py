import argparse
import os
import sys
import tqdm
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from obtain_annotations import main as obtain_annotations_main


def main(main_proteomes_dir, output_dir, go_file, proteomes_ids_file):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'annotations'), exist_ok=True)
    with open(proteomes_ids_file, 'r') as f:
        for line in tqdm.tqdm(f):
            if line.startswith('TaxonID'):
                continue
            parts = line.strip().split('\t')
            taxon_id = parts[0]
            proteome_id = parts[1]
            
            # Search under Archaea, Bacteria, Eukaryota, and Viruses directories
            proteome_subdir = None
            for domain in ['Bacteria', 'Eukaryota', 'Archaea', 'Viruses']:
                potential_path = os.path.join(main_proteomes_dir, domain, proteome_id)
                if os.path.exists(potential_path):
                    proteome_subdir = potential_path
                    break
            
            if proteome_subdir is None:
                print(f"Warning: Proteome {proteome_id} not found in any domain directory", file=sys.stderr)
                continue
            
            proteome_file = os.path.join(proteome_subdir, proteome_id + '_' + taxon_id + '.dat.gz')
            if not os.path.exists(proteome_file):
                print(f"Warning: Proteome file {proteome_file} does not exist", file=sys.stderr)
                continue
            output_file = os.path.join(output_dir, 'annots_taxon_' + taxon_id + '.tsv')
            obtain_annotations_main(proteome_file, output_file, go_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--main_proteomes_dir', type=str, default="${DATA_DIR}/uniprot_reference_proteomes")
    parser.add_argument('--output_dir', type=str, default="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic")
    parser.add_argument('--go_file', type=str, default="data/go-basic.obo")
    parser.add_argument('--proteomes_ids_file', type=str, default="data/uniprot_proteomes_ids.tsv")
    args = parser.parse_args()

    main(args.main_proteomes_dir, args.output_dir, args.go_file, args.proteomes_ids_file)