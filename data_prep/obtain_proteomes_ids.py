"""
Given swissprot_exp_2025_08_proteomes.pkl and viruses.txt, obtain the proteomes ids excluding the viruses.
"""

import pandas as pd
import numpy as np
import os
import sys
import argparse

def load_viruses_file(viruses_file):
    viruses = set()
    with open(viruses_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not line.startswith('#'):
                viruses.add(parts[0].strip())
    return list(viruses)


def load_available_reference_proteomes(uniprot_proteomes_readme_file):
    proteomes_to_tax_id_and_entries_count_and_domain = dict()
    record = False
    with open(uniprot_proteomes_readme_file, 'r') as f:
        for line in f:
            if line.startswith('Proteome_ID	Tax_ID	OSCODE	SUPERREGNUM	#(1)	#(2)	#(3)	Species Name'):
                record = True
            elif record:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    record = False
                    break
                proteome_id = parts[0].strip()
                tax_id = parts[1].strip()
                domain = parts[3].strip()
                entries_count = int(parts[4].strip())
                proteomes_to_tax_id_and_entries_count_and_domain[proteome_id] = (tax_id, entries_count, domain)
    return proteomes_to_tax_id_and_entries_count_and_domain


def obtain_proteomes_ids(swissprot_file, viruses_file, uniprot_proteomes_readme_file):
    viruses = load_viruses_file(viruses_file)
    reference_proteomes_to_tax_id_and_entries_count_and_domain = load_available_reference_proteomes(uniprot_proteomes_readme_file)
    swissprot = pd.read_pickle(swissprot_file)
    filtered_swissprot = swissprot[~swissprot['orgs'].isin(viruses)]
    number_of_lines_without_viruses = filtered_swissprot.shape[0]

    filtered_swissprot = filtered_swissprot[filtered_swissprot['proteomes'].apply(lambda x: len(x) > 0 and x is not None)]
    print(f"dropped {number_of_lines_without_viruses - filtered_swissprot.shape[0]} lines without proteomes")
    print(f"Number of lines in the filtered swissprot: {filtered_swissprot.shape[0]}")

    # proteomes ids are stored in 'proteomes' column as an array of strings
    proteomes_ids = filtered_swissprot['proteomes'].explode().dropna().unique()

    # want to make sure that proteomes of same organisms are common
    unique_organisms = set()
    common_organisms = dict()
    organism_to_proteomes = dict()
    for organism in filtered_swissprot['orgs'].unique():
        organism_proteomes = filtered_swissprot[filtered_swissprot['orgs'] == organism]['proteomes'].explode().dropna().unique()
        if len(organism_proteomes) > 1:
            common_organisms[organism] = organism_proteomes
        
        if len(organism_proteomes) >= 1:
            unique_organisms.add(organism)
            organism_to_proteomes[organism] = organism_proteomes
    
    print(f"Unique organisms with proteomes: {len(unique_organisms)}")
    print(f"Organisms with multiple proteomes: {len(common_organisms)}")

    # get count of lines per proteome for the common organisms (order by organism id and print each organism and its proteomes count)
    common_organisms_proteomes_counts = filtered_swissprot[filtered_swissprot['orgs'].isin(common_organisms)]['proteomes'].explode().dropna().value_counts()
    i = 0
    for organism, proteomes in common_organisms.items():
        i += 1
        for proteome in proteomes:
            print(f"{i:02d}. Organism {organism}: {proteome} {common_organisms_proteomes_counts[proteome]}")
        print("------")

    # pick reference proteome for each organism
    unique_organism_proteome = dict()
    proteome_counts = filtered_swissprot['proteomes'].explode().dropna().value_counts()
    for organism, proteomes in organism_to_proteomes.items():
        found_reference_proteome = False
        ordered_proteomes = sorted(proteomes, key=lambda x: proteome_counts[x], reverse=True)
        for proteome in ordered_proteomes:
            if proteome in reference_proteomes_to_tax_id_and_entries_count_and_domain:
                tax_id, entries_count, domain = reference_proteomes_to_tax_id_and_entries_count_and_domain[proteome]
                found_reference_proteome = True
                break
        if not found_reference_proteome:
            print(f"Warning: No reference proteome found for organism {organism}")
            continue
        unique_organism_proteome[organism] = (proteome, proteome_counts[proteome], entries_count, domain)

    print(f"Number of unique organism proteomes: {len(unique_organism_proteome)}")

    return unique_organism_proteome

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--swissprot_file', type=str, default='data/swissprot_exp_2025_10_proteomes.pkl')
    parser.add_argument('--viruses_file', type=str, default='viruses.txt')
    parser.add_argument('--output_file', type=str, default='data/uniprot_proteomes_ids.tsv')
    parser.add_argument('--uniprot_proteomes_readme_file', type=str, default='data/uniprot_proteomes_readme.txt')
    args = parser.parse_args()

    proteomes_ids = obtain_proteomes_ids(args.swissprot_file, args.viruses_file, args.uniprot_proteomes_readme_file)
    
    with open(args.output_file, 'w') as f:
        f.write(f"TaxonID\tProteomeID\tEXPCount\tEntriesCount\tDomain\n")
        # sort by counts descending
        proteomes_ids_sorted = sorted(proteomes_ids.items(), key=lambda x: x[1][1], reverse=True)
        for organism, (proteome_id, counts, entries_count, domain) in proteomes_ids_sorted:
            f.write(f"{organism}\t{proteome_id}\t{counts}\t{entries_count}\t{domain}\n")
    print(f"Successfully saved {len(proteomes_ids_sorted)} organism proteomes to {args.output_file}")