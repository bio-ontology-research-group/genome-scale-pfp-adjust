#!/usr/bin/env python
"""
Extract proteins from UniProt/SwissProt .dat file based on integration date.
This script parses the SwissProt .dat file and extracts protein names for proteins
that were integrated into UniProtKB/Swiss-Prot after a specified date.
"""

import os
import sys
import argparse
import gzip
import logging
from datetime import datetime
from typing import List, Tuple, Optional, Set

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)


def parse_uniprot_date(date_str: str) -> Optional[datetime]:
    """
    Convert UniProt date format (DD-MMM-YYYY) to datetime object.
    
    Args:
        date_str: Date string in format like "23-MAY-2024" or "02-NOV-2023"
    
    Returns:
        datetime object or None if parsing fails
    """
    try:
        # UniProt uses format like "23-MAY-2024"
        return datetime.strptime(date_str.strip(), "%d-%b-%Y")
    except ValueError as e:
        logging.warning(f"Failed to parse date '{date_str}': {e}")
        return None


def load_proteome_taxon_ids(proteomes_file: str) -> Set[str]:
    """
    Load taxon IDs from the UniProt proteomes TSV file.
    
    Args:
        proteomes_file: Path to uniprot_proteomes_ids.tsv file
    
    Returns:
        Set of taxon IDs as strings
    """
    taxon_ids = set()
    
    if not os.path.exists(proteomes_file):
        logging.warning(f"Proteomes file not found: {proteomes_file}")
        return taxon_ids
    
    logging.info(f"Loading proteome taxon IDs from: {proteomes_file}")
    
    with open(proteomes_file, 'r') as f:
        # Skip header
        header = f.readline()
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if parts:
                taxon_id = parts[0].strip()
                taxon_ids.add(taxon_id)
    
    logging.info(f"Loaded {len(taxon_ids)} taxon IDs from proteomes file")
    return taxon_ids


def parse_dat_file(file_path: str, cutoff_date: datetime, proteome_taxon_ids: Set[str]) -> List[Tuple[str, str, str, str, bool]]:
    """
    Parse UniProt/SwissProt .dat file and extract proteins integrated after cutoff date.
    
    Args:
        file_path: Path to SwissProt .dat file (can be gzipped)
        cutoff_date: Only return proteins integrated after this date
        proteome_taxon_ids: Set of taxon IDs from proteomes file for matching
    
    Returns:
        List of tuples (protein_id, protein_name, integration_date_str, matched_taxon_ids, is_matched)
    """
    results = []
    
    # Open file (handle both gzipped and plain text)
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    logging.info(f"Parsing SwissProt file: {file_path}")
    logging.info(f"Looking for proteins integrated after: {cutoff_date.strftime('%d-%b-%Y')}")
    
    with open_func(file_path, mode) as f:
        prot_id = ''
        prot_name = ''
        integration_date_str = ''
        integration_date = None
        taxon_ids = []
        
        line_count = 0
        protein_count = 0
        
        for line in f:
            line_count += 1
            if line_count % 10000000 == 0:
                logging.info(f"Processed {line_count} lines, found {len(results)} matching proteins")
            
            # Split line by multiple spaces (UniProt format uses 3+ spaces as delimiter)
            items = line.strip().split('   ')
            
            if items[0] == 'ID' and len(items) > 1:
                # Save previous protein if it matches criteria
                if prot_id != '' and integration_date and integration_date > cutoff_date:
                    # Match taxon IDs against proteomes
                    matched_taxons = [tid for tid in taxon_ids if tid in proteome_taxon_ids]
                    matched_taxon_str = ';'.join(matched_taxons) if matched_taxons else 'NA'
                    is_matched = len(matched_taxons) > 0
                    results.append((prot_id, prot_name, integration_date_str, matched_taxon_str, is_matched))
                
                # Start new protein
                prot_id = items[1].split()[0]  # Get just the ID part
                prot_name = ''
                integration_date_str = ''
                integration_date = None
                taxon_ids = []
                protein_count += 1
            
            elif items[0] == 'AC' and len(items) > 1:
                # Could use accession if needed, but we're using ID for now
                pass
            
            elif items[0] == 'OX' and len(items) > 1:
                # Parse OX line to extract taxon IDs
                # Format: "NCBI_TaxID=654924;"
                ox_line = items[1]
                if 'NCBI_TaxID=' in ox_line:
                    # Extract the taxon ID(s) - can have multiple separated by commas
                    taxid_part = ox_line.split('NCBI_TaxID=')[1].split(';')[0]
                    # Handle potential multiple IDs (though typically one)
                    for tid in taxid_part.split(','):
                        tid = tid.strip()
                        if tid:
                            taxon_ids.append(tid)
            
            elif items[0] == 'DT' and len(items) > 1:
                # Parse date lines
                dt_line = items[1]
                if 'integrated into UniProtKB/Swiss-Prot' in dt_line:
                    # Extract date from line like "23-MAY-2024, integrated into UniProtKB/Swiss-Prot."
                    date_part = dt_line.split(',')[0].strip()
                    integration_date_str = date_part
                    integration_date = parse_uniprot_date(date_part)
            
            elif items[0] == 'DE' and len(items) > 1:
                # Parse protein names - prioritize RecName
                de_line = items[1]
                if de_line.startswith('RecName: Full='):
                    name = de_line.split('Full=')[1].split(';')[0].strip()
                    prot_name = name
                elif de_line.startswith('AltName: Full=') and not prot_name:
                    # Use AltName only if RecName not found
                    name = de_line.split('Full=')[1].split(';')[0].strip()
                    prot_name = name
        
        # Save last protein if it matches criteria
        if prot_id != '' and integration_date and integration_date > cutoff_date:
            # Match taxon IDs against proteomes
            matched_taxons = [tid for tid in taxon_ids if tid in proteome_taxon_ids]
            matched_taxon_str = ';'.join(matched_taxons) if matched_taxons else 'NA'
            is_matched = len(matched_taxons) > 0
            results.append((prot_id, prot_name, integration_date_str, matched_taxon_str, is_matched))
    
    logging.info(f"Parsed {protein_count} total proteins")
    logging.info(f"Found {len(results)} proteins integrated after {cutoff_date.strftime('%d-%b-%Y')}")
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract proteins from SwissProt .dat file by integration date'
    )
    parser.add_argument(
        '--dat-file',
        default='data/uniprot_sprot_2025_10.dat',
        help='Path to UniProt/SwissProt .dat file (can be gzipped)'
    )
    parser.add_argument(
        '--cutoff-date',
        required=True,
        help='Cutoff date in format DD-MMM-YYYY (e.g., "23-MAY-2024" or "02-NOV-2023")'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output TSV file path'
    )
    parser.add_argument(
        '--proteomes-file',
        default='data/uniprot_proteomes_ids.tsv',
        help='Path to UniProt proteomes TSV file with taxon IDs'
    )
    
    args = parser.parse_args()

    if args.output is None:
        args.output = f"proteins_by_date_{args.cutoff_date}.tsv"
    
    # Validate and parse cutoff date
    cutoff_date = parse_uniprot_date(args.cutoff_date)
    if not cutoff_date:
        logging.error(f"Invalid cutoff date format: {args.cutoff_date}")
        logging.error("Expected format: DD-MMM-YYYY (e.g., '23-MAY-2024')")
        sys.exit(1)
    
    # Check if input file exists
    if not os.path.exists(args.dat_file):
        logging.error(f"Input file not found: {args.dat_file}")
        sys.exit(1)
    
    # Load proteome taxon IDs
    proteome_taxon_ids = load_proteome_taxon_ids(args.proteomes_file)
    
    # Parse the file
    results = parse_dat_file(args.dat_file, cutoff_date, proteome_taxon_ids)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")
    
    # Write results to TSV
    with open(args.output, 'w') as f:
        # Write header
        f.write("protein_id\tprotein_name\tintegration_date\ttaxon_id\ttaxon_id_matched\n")
        
        # Write data
        for prot_id, prot_name, integration_date, taxon_id, is_matched in results:
            # Handle missing protein names
            if not prot_name:
                prot_name = "N/A"
            # Convert boolean to string
            matched_str = "True" if is_matched else "False"
            f.write(f"{prot_id}\t{prot_name}\t{integration_date}\t{taxon_id}\t{matched_str}\n")
    
    logging.info(f"Results saved to: {args.output}")
    logging.info(f"Total proteins written: {len(results)}")


if __name__ == '__main__':
    main()
