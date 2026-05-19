"""
Run the adjustment pipeline.

Steps:

1. Load constraints, go hierarchy, taxon hierarchy, ncbitaxon hierarchy.
2. Run taxon consistency adjustment.
3. Run complex coherence adjustment.
"""

from typing import Dict
import os
import argparse
import json
import sys
from tqdm import tqdm
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# from taxon_consistency.adjust_ortools import main as taxon_consistency_adjust
from taxon_consistency.adjust_ortools import load_constraints, load_go_hierarchy, load_taxon_hierarchy, adjust_per_taxon

from complex_coherence.adjust_ortools import main as complex_coherence_adjust



def get_taxons_and_folds(predictions_dir: str):
    # Get all taxons and folds from the predictions directory
    taxons_str = []
    folds_str = []
    for file in os.listdir(predictions_dir):
        if file.startswith('predictions_fold_') and '_taxon_' in file and file.endswith('.tsv'):
            taxons_str.append(file.split('_taxon_')[1].split('.')[0].split('_')[0])
            folds_str.append(file.split('fold_')[1].split('_')[0])
    return taxons_str, folds_str

def verify_predictions_format(predictions_file: str):
    # Verify predictions are in the correct format 
    # predictions_file is a tsv file with the following columns:
    # 1. protein_id
    # 2. GO:XXXXXX|score
    # 3. GO:YYYYY|score
    # 4. ...
    # The protein_id is the first column and the GO terms are the remaining columns.
    # The GO terms are separated by '|' and the score is the second part.
    # The score is a float number.
    for line in open(predictions_file, 'r'):
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            return False
        for part in parts[1:]:
            if '|' not in part:
                return False
    return True



def run_adjustment_pipeline(
    predictions_dir: str,
    optimized_dir: str,
    fold_thresholds: Dict[str, float],
    constraints_file: str,
    go_hierarchy_file: str,
    taxon_hierarchy_file: str,
    ncbitaxon_hierarchy_file: str,
    complexes_file: str,
    optimized: bool = False,
    complex_coherence: bool = False,
    top_k: int = None,
    start_index: int = 0,
    end_index: int = None,
    skip_existing: bool = False,
    provide_taxon_id: bool = False,
    taxon_ids_file: str = None,
):
    """
    Run the adjustment pipeline.
    """
    in_taxon_constraints, never_in_taxon_constraints = load_constraints(constraints_file)
    go_hierarchy = load_go_hierarchy(go_hierarchy_file)
    taxon_is_a_hierarchy, taxon_disjoint_from_hierarchy, taxon_union_of_hierarchy = load_taxon_hierarchy(taxon_hierarchy_file, ncbitaxon_hierarchy_file)
    print(f"Loaded {len(in_taxon_constraints)} 'only_in_taxon' constraints")
    print(f"Loaded {len(never_in_taxon_constraints)} 'never_in_taxon' constraints")
    print(f"Loaded {len(go_hierarchy)} GO hierarchy relationships")
    print(f"Loaded {len(taxon_is_a_hierarchy)} taxon is_a hierarchy relationships")
    print(f"Loaded {len(taxon_disjoint_from_hierarchy)} taxon disjoint_from hierarchy relationships")
    print(f"Loaded {len(taxon_union_of_hierarchy)} taxon union_of hierarchy relationships")
    print("--------------------------------")
    taxons_str, folds_str = get_taxons_and_folds(predictions_dir)

    # Sort for consistent ordering across parallel runs
    taxon_fold_pairs = sorted(zip(taxons_str, folds_str))
    taxons_str = [t for t, f in taxon_fold_pairs]
    folds_str = [f for t, f in taxon_fold_pairs]

    # Filter by taxon IDs file if provided
    if taxon_ids_file is not None:
        with open(taxon_ids_file, 'r') as f:
            allowed_taxon_ids = {line.strip() for line in f if line.strip()}
        taxon_fold_pairs = [(t, f) for t, f in taxon_fold_pairs if t in allowed_taxon_ids]
        taxons_str = [t for t, f in taxon_fold_pairs]
        folds_str = [f for t, f in taxon_fold_pairs]

    total_taxons = len(taxons_str)

    # Apply index range for parallel processing
    if end_index is None:
        end_index = len(taxons_str)
    taxons_str = taxons_str[start_index:end_index]
    folds_str = folds_str[start_index:end_index]

    print(f"Processing taxons {start_index} to {end_index} (out of {total_taxons} total)")

    num_flips = 0
    num_annotated_predictions = 0
    complex_coherence_num_flips = 0
    complex_coherence_num_annotated_predictions = 0
    for i, taxon_str in tqdm(enumerate(taxons_str)):
        fold_str = folds_str[i]

        genome_predictions_file = os.path.join(predictions_dir, f'predictions_fold_{fold_str}_taxon_{taxon_str}.tsv')

        genome_optimized_output_file = os.path.join(optimized_dir, f'optimized_fold_{fold_str}_taxon_{taxon_str}.tsv')
        genome_threshold = fold_thresholds[fold_str]

        # Skip if output already exists and skip_existing is True
        if skip_existing and os.path.exists(genome_optimized_output_file):
            print(f"Skipping {genome_optimized_output_file} (already exists)")
            continue

        # Run adjustments scripts
        # Taxon consistency adjustment

        # taxon_consistency_adjust(
        #     predictions_file=genome_predictions_file,
        #     constraints_file=constraints_file,
        #     go_hierarchy_file=go_hierarchy_file,
        #     taxon_hierarchy_file=taxon_hierarchy_file,
        #     ncbitaxon_hierarchy_file=ncbitaxon_hierarchy_file,
        #     output_file=genome_optimized_output_file,
        #     threshold=genome_threshold,
        # )

        # Taxon consistency adjustment (faster to load files once)
        genome_total_flips, genome_num_annotated_predictions = adjust_per_taxon(
            predictions_file=genome_predictions_file,
            in_taxon_constraints=in_taxon_constraints,
            never_in_taxon_constraints=never_in_taxon_constraints,
            go_hierarchy=go_hierarchy,
            taxon_is_a_hierarchy=taxon_is_a_hierarchy,
            taxon_disjoint_from_hierarchy=taxon_disjoint_from_hierarchy,
            taxon_union_of_hierarchy=taxon_union_of_hierarchy,
            output_file=genome_optimized_output_file,
            threshold=genome_threshold,
            taxon_id=taxon_str if provide_taxon_id else None,
        )
        num_flips += genome_total_flips
        num_annotated_predictions += genome_num_annotated_predictions



        # Complex coherence adjustment
        if complex_coherence:
            genome_complex_coherence_total_flips, genome_complex_coherence_num_annotated_predictions = complex_coherence_adjust(
                predictions_file=genome_optimized_output_file,
                complexes_file=complexes_file,
                go_hierarchy_file=go_hierarchy_file,
                output_file=genome_optimized_output_file,
                threshold=genome_threshold,
                optimized=optimized,
                top_k=top_k,
            )

            complex_coherence_num_flips += genome_complex_coherence_total_flips
            complex_coherence_num_annotated_predictions += genome_complex_coherence_num_annotated_predictions


    # print summary of the results
    print("--------------------------------")
    print("Summary of the results:")
    print("--------------------------------")
    print(f"Total taxons: {len(taxons_str)}")
    print(f"Total folds: {len(folds_str)}")
    print(f"Total flips (taxon consistency): {num_flips}")
    print(f"Total annotated predictions (taxon consistency): {num_annotated_predictions}")
    print(f"Percentage of flipped predictions (taxon consistency): {num_flips / num_annotated_predictions * 100:.2f}%")
    if complex_coherence:
        print(f"Total flips (complex coherence): {complex_coherence_num_flips}")
        print(f"Total annotated predictions (complex coherence): {complex_coherence_num_annotated_predictions}")
        print(f"Percentage of flipped predictions (complex coherence): {complex_coherence_num_flips / complex_coherence_num_annotated_predictions * 100:.2f}%")
    print("--------------------------------")
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run the adjustment pipeline.')
    parser.add_argument('--predictions_dir', required=True, help='Path to predictions directory')
    parser.add_argument('--optimized_dir', default=None, help='Path to optimized directory')
    parser.add_argument('--fold_thresholds', default=None, help='Path to fold thresholds file')
    parser.add_argument('--constraints_file', default="data/go_taxon_constraints_extracted_obo.tsv", help='Path to constraints file')
    parser.add_argument('--go_hierarchy_file', default="data/go_hierarchy.tsv", help='Path to go hierarchy file')
    parser.add_argument('--taxon_hierarchy_file', default="data/taxon_hierarchy.tsv", help='Path to taxon hierarchy file')
    parser.add_argument('--ncbitaxon_hierarchy_file', default="data/ncbitaxon_hierarchy.tsv", help='Path to ncbitaxon hierarchy file')
    
    parser.add_argument('--complexes_file', default="complex_coherence/protein_complexes.tsv", help='Path to complexes file')
    parser.add_argument('--optimized', action='store_true', help='Use optimized constraints')
    parser.add_argument('--complex_coherence', action='store_true', help='Use complex coherence adjustment')
    parser.add_argument('--top_k', type=int, default=None, help='Top-k proteins to participate in each complex (default: None)')
    parser.add_argument('--start_index', type=int, default=0, help='Starting index in taxon list (inclusive, for parallel processing)')
    parser.add_argument('--end_index', type=int, default=None, help='Ending index in taxon list (exclusive, for parallel processing)')
    parser.add_argument('--skip_existing', action='store_true', help='Skip taxons whose output file already exists')
    parser.add_argument('--provide_taxon_id', action='store_true', help='Provide taxon_id to the solver (default: None)')
    parser.add_argument('--taxon_ids_file', default=None, help='Path to txt file with taxon IDs to process (one per line). Only taxons present in both this file and predictions_dir will be processed.')
    parser.add_argument('--default_threshold', type=float, default=0.5, help='Default threshold for fold thresholds')
    args = parser.parse_args()

    if args.optimized_dir is None:
        args.optimized_dir = args.predictions_dir.replace('predictions', 'optimized')
    os.makedirs(args.optimized_dir, exist_ok=True)

    if args.fold_thresholds is not None:
        with open(args.fold_thresholds, 'r') as f:
            args.fold_thresholds = json.load(f)
    else:
        # try to load fold_thresholds.json from the parent of predictions directory
        fold_thresholds_file = os.path.join(os.path.dirname(args.predictions_dir), 'fold_thresholds.json')
        if os.path.exists(fold_thresholds_file):
            with open(fold_thresholds_file, 'r') as f:
                args.fold_thresholds = json.load(f)
        else:
            # set default fold thresholds to 0.5
            folds = [int(fold.split('_')[2]) for fold in os.listdir(args.predictions_dir) if fold.startswith('predictions_fold_') and fold.endswith('.tsv')]
            args.fold_thresholds = {f"{fold:02d}": args.default_threshold for fold in folds}
        
    print(f"Predictions directory: {args.predictions_dir}")
    print(f"Optimized directory: {args.optimized_dir}")
    print(f"Fold thresholds: {args.fold_thresholds}")
    
    
    run_adjustment_pipeline(
        predictions_dir=args.predictions_dir,
        optimized_dir=args.optimized_dir,
        fold_thresholds=args.fold_thresholds,
        constraints_file=args.constraints_file,
        go_hierarchy_file=args.go_hierarchy_file,
        taxon_hierarchy_file=args.taxon_hierarchy_file,
        ncbitaxon_hierarchy_file=args.ncbitaxon_hierarchy_file,
        complexes_file=args.complexes_file,
        optimized=args.optimized,
        complex_coherence=args.complex_coherence,
        top_k=args.top_k,
        start_index=args.start_index,
        end_index=args.end_index,
        skip_existing=args.skip_existing,
        provide_taxon_id=args.provide_taxon_id,
        taxon_ids_file=args.taxon_ids_file,
    )