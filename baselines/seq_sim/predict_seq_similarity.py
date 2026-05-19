"""
Sequence similarity based functional annotator.

Predicts GO functions for all proteins of test organisms (test_organisms.txt) using
Diamond BLAST against a database built from annotated train organism proteins
(per-taxon annotation TSV files + UniProt reference proteome FASTAs).

This is a one-shot (non-CV) counterpart to seq_similarity.py: train/test splits are
provided as explicit organism lists rather than cross-validation folds.

Usage example:
    python predict_seq_similarity.py \\
        --subontology mf \\
        --train-organisms-file splits/heldout/train_organisms.txt \\
        --test-organisms-file  splits/heldout/test_organisms.txt \\
        --annotations-dir      ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \\
        --main-proteomes-dir   ${DATA_DIR}/uniprot_reference_proteomes \\
        --proteomes-ids-file   data/uniprot_proteomes_ids.tsv \\
        --output-dir           seq_similarity_heldout_results_mf
"""

import os
import json
import gc
import sys
from typing import Dict, List, Set

import click as ck
import numpy as np
import pandas as pd

try:
    import psutil
except ImportError:
    psutil = None

# Allow imports from the parent directory (deepgo modules)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', 'GAEF')))

# Import shared utilities from seq_similarity (safe: main() is guarded by __name__)
import seq_similarity as ss


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_organisms_from_file(file_path: str) -> Set[str]:
    """
    Load organism taxids from a plain-text file (one taxid per line).
    Lines that are blank or start with '#' are skipped.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Organisms file not found: {file_path}")

    organisms: Set[str] = set()
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            taxid = line.strip()
            if not taxid or taxid.startswith('#'):
                continue
            if not taxid.isdigit():
                ck.echo(f"[WARNING] Invalid taxid '{taxid}' on line {line_num} in {file_path}")
                continue
            organisms.add(taxid)

    ck.echo(f"[INFO] Loaded {len(organisms)} organisms from {file_path}")
    return organisms


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

ALL_SUBONTOLOGIES = ['mf', 'bp', 'cc']


def _resolve_terms_file(go_obo: str, subont: str) -> str:
    """Locate terms.pkl for *subont*, trying beside the OBO file then beside the script."""
    go_obo_dir = os.path.dirname(os.path.abspath(go_obo))
    candidate = os.path.join(go_obo_dir, subont, 'terms.pkl')
    if os.path.exists(candidate):
        return candidate
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    fallback = os.path.join(repo_root, 'data', subont, 'terms.pkl')
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(
        f"[ERROR] terms.pkl not found for sub-ontology '{subont}'. "
        f"Tried:\n  {candidate}\n  {fallback}"
    )


@ck.command()
@ck.option('--subontology', default='cc',
           type=ck.Choice(['mf', 'bp', 'cc', 'all']),
           help="GO sub-ontology to predict. Use 'all' to run mf+bp+cc with a single BLAST.")
@ck.option('--go-obo', default='data/go-basic.obo',
           help='GO OBO file (used to locate terms.pkl).')
@ck.option('--train-organisms-file',
           default='splits/heldout/train_organisms.txt',
           help='File listing train organism taxids (one per line).')
@ck.option('--test-organisms-file',
           default='splits/heldout/test_organisms.txt',
           help='File listing test organism taxids (one per line).')
@ck.option('--annotations-dir',
           default='${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic',
           help='Directory containing per-taxon annotation TSV files (annots_taxon_<id>.tsv).')
@ck.option('--main-proteomes-dir',
           default='${DATA_DIR}/uniprot_reference_proteomes',
           help='UniProt reference proteomes base directory.')
@ck.option('--proteomes-ids-file',
           default='data/uniprot_proteomes_ids.tsv',
           help='Proteomes IDs mapping file (TaxonID <TAB> ProteomeID).')
@ck.option('--output-dir', 
           default='${DATA_DIR}/swissprot_proteomes_folds/seq_sim_heldout_results',
           help='Output directory.')
@ck.option('--evalue', type=float, default=1e-5,
           help='Diamond BLAST E-value threshold.')
@ck.option('--max-target-seqs', type=int, default=50,
           help='Maximum target sequences for Diamond BLAST.')
@ck.option('--threads', type=int, default=None,
           help='Number of threads for Diamond BLAST (default: auto-detect).')
@ck.option('--memory-limit', default=None,
           help='Memory limit for Diamond BLAST (e.g., 32G, 64G).')
@ck.option('--max-workers', type=int, default=None,
           help='Maximum parallel workers for FASTA processing (default: auto).')
@ck.option('--chunk-size', type=int, default=1_000_000,
           help='Chunk size for parsing Diamond results.')
@ck.option('--batch-size', type=int, default=500_000,
           help='Batch size for splitting large query FASTA files.')
@ck.option('--sparsity-threshold', type=float, default=0.001,
           help='Minimum score to include in output TSVs.')
def main(
    subontology, go_obo,
    train_organisms_file, test_organisms_file, annotations_dir,
    main_proteomes_dir, proteomes_ids_file,
    output_dir,
    evalue, max_target_seqs, threads, memory_limit,
    max_workers, chunk_size, batch_size, sparsity_threshold,
):
    """Sequence-similarity annotator: train organisms → Diamond DB → predict test organisms.

    When --subontology all is used, Diamond BLAST runs only once and GO term
    assignments are computed separately per sub-ontology (mf, bp, cc).
    Outputs land in predictions/test_mf/, predictions/test_bp/, predictions/test_cc/.
    """

    # ------------------------------------------------------------------
    # 0. Validate / normalise paths
    # ------------------------------------------------------------------
    main_proteomes_dir   = os.path.abspath(main_proteomes_dir)
    proteomes_ids_file   = os.path.abspath(proteomes_ids_file)
    train_organisms_file = os.path.abspath(train_organisms_file)
    test_organisms_file  = os.path.abspath(test_organisms_file)
    annotations_dir      = os.path.abspath(annotations_dir)
    output_dir           = os.path.abspath(output_dir)

    for path, label in [
        (main_proteomes_dir,   'main_proteomes_dir'),
        (proteomes_ids_file,   'proteomes_ids_file'),
        (train_organisms_file, 'train_organisms_file'),
        (test_organisms_file,  'test_organisms_file'),
        (annotations_dir,      'annotations_dir'),
    ]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"[ERROR] {label} not found: {path}")

    os.makedirs(output_dir, exist_ok=True)

    # Resolve which sub-ontologies to run (fail fast on missing terms files)
    subontologies: List[str] = ALL_SUBONTOLOGIES if subontology == 'all' else [subontology]
    terms_files: Dict[str, str] = {}
    for subont in subontologies:
        terms_files[subont] = _resolve_terms_file(go_obo, subont)
        ck.echo(f"[INFO] terms.pkl for '{subont}': {terms_files[subont]}")

    # ------------------------------------------------------------------
    # 1. Load organism sets
    # ------------------------------------------------------------------
    train_organisms = load_organisms_from_file(train_organisms_file)
    test_organisms  = load_organisms_from_file(test_organisms_file)

    overlap = train_organisms & test_organisms
    if overlap:
        ck.echo(f"[WARNING] {len(overlap)} organisms appear in both train and test sets: "
                f"{sorted(overlap)}")

    ck.echo(f"[INFO] Train organisms: {len(train_organisms)}, "
            f"Test organisms: {len(test_organisms)}")

    # ------------------------------------------------------------------
    # 2. Load proteome mapping
    # ------------------------------------------------------------------
    proteome_mapping = ss.load_proteome_mapping(proteomes_ids_file)

    # ------------------------------------------------------------------
    # 3. Load train annotations + sequences from per-taxon files
    #    ss.create_dataframe_from_annotations_and_sequences reads
    #    annots_taxon_<id>.tsv from annotations_dir and the matching
    #    .fasta.gz from main_proteomes_dir for every train organism.
    #    Annotations cover ALL GO terms; sub-ontology filtering happens later.
    # ------------------------------------------------------------------
    ck.echo("[INFO] Loading train annotations and sequences...")
    train_df = ss.create_dataframe_from_annotations_and_sequences(
        annotations_dir, main_proteomes_dir, proteome_mapping, train_organisms
    )
    ck.echo(f"[INFO] Loaded {len(train_df)} train proteins")

    if train_df.empty:
        raise ValueError(
            "[ERROR] No train proteins loaded. "
            "Check --annotations-dir and --train-organisms-file."
        )

    # Keep only annotated proteins for the DB and the annotations dict
    train_df = train_df[train_df['prop_annotations'].apply(len) > 0].copy()
    ck.echo(f"[INFO] {len(train_df)} train proteins have GO annotations")

    train_annotations: Dict[str, Set[str]] = {
        row.proteins: set(row.prop_annotations)
        for row in train_df.itertuples()
    }
    n_train_annotations = len(train_annotations)

    train_df_with_seq = train_df[train_df['sequences'] != ''].copy()
    ck.echo(f"[INFO] {len(train_df_with_seq)}/{len(train_df)} annotated train proteins have sequences")

    if train_df_with_seq.empty:
        raise ValueError("[ERROR] No train protein sequences found. Check --main-proteomes-dir.")

    del train_df   # keep only the sequence-filtered copy going forward
    gc.collect()

    # ------------------------------------------------------------------
    # 4. Create Diamond database  (sub-ontology agnostic)
    # ------------------------------------------------------------------
    ck.echo("[INFO] Creating Diamond database from train proteins...")
    db_file = ss.create_diamond_database(
        train_df_with_seq,
        pd.DataFrame(),   # no supplemental model-organism sequences
        output_dir,
        fold_id=1,
    )
    ck.echo(f"[INFO] Diamond database: {db_file}")

    del train_df_with_seq
    gc.collect()

    # ------------------------------------------------------------------
    # 5. Create query FASTA for test organisms  (sub-ontology agnostic)
    # ------------------------------------------------------------------
    ck.echo("[INFO] Creating query FASTA for test organisms...")
    query_fasta = ss.create_query_fasta(
        test_organisms, main_proteomes_dir, output_dir,
        fold_id=1, proteome_mapping=proteome_mapping,
        max_workers=max_workers,
    )
    ck.echo(f"[INFO] Query FASTA: {query_fasta}")

    # ------------------------------------------------------------------
    # 6. Run Diamond BLAST  (sub-ontology agnostic — done exactly once)
    # ------------------------------------------------------------------
    diamond_output = os.path.join(output_dir, 'diamond_fold_01.txt')
    ck.echo("[INFO] Running Diamond BLAST...")
    success = ss.run_diamond_blast(
        query_fasta, db_file, diamond_output,
        evalue, max_target_seqs, threads, memory_limit, batch_size,
    )

    if not success:
        raise RuntimeError("[ERROR] Diamond BLAST failed. Check logs above for details.")

    ck.echo(f"[INFO] Diamond output: {diamond_output}")

    # ------------------------------------------------------------------
    # 7. Parse Diamond results  (sub-ontology agnostic — done exactly once)
    # ------------------------------------------------------------------
    ck.echo("[INFO] Parsing Diamond results...")
    diamond_results = ss.parse_diamond_results(diamond_output, chunk_size)
    ck.echo(f"[INFO] Diamond results: {len(diamond_results)} query proteins with hits")

    # ------------------------------------------------------------------
    # 8. Build protein → taxon mapping  (sub-ontology agnostic)
    # ------------------------------------------------------------------
    ck.echo("[INFO] Building protein → taxon mapping for test organisms...")
    protein_to_taxid = ss.get_protein_to_taxid_mapping(
        main_proteomes_dir, test_organisms, proteome_mapping
    )

    # ------------------------------------------------------------------
    # 9–11. Per-sub-ontology: compute predictions → organise → write TSVs
    # ------------------------------------------------------------------
    all_summaries: List[Dict] = []

    for subont in subontologies:
        ck.echo(f"\n[INFO] === Processing sub-ontology: {subont.upper()} ===")

        # Load terms for this sub-ontology
        terms_array = pd.read_pickle(terms_files[subont])['gos'].values.flatten()
        all_terms: List[str] = list(terms_array)
        ck.echo(f"[INFO] Loaded {len(all_terms)} GO terms for '{subont}'")

        # Compute bitscore-weighted GO scores (filtered to this ontology's terms)
        ck.echo(f"[INFO] Computing GO predictions for '{subont}'...")
        predictions = ss.compute_diamond_predictions(
            diamond_results, train_annotations, all_terms
        )
        ck.echo(f"[INFO] Predictions generated for {len(predictions)} proteins")

        # Organise by taxon and write TSVs
        predictions_by_taxon = ss.organize_predictions_by_taxon(
            predictions, protein_to_taxid, test_organisms
        )

        ck.echo(f"[INFO] Writing per-taxon TSVs for '{subont}' "
                f"(sparsity_threshold={sparsity_threshold})...")
        written_files = ss.write_taxon_tsvs(
            predictions_by_taxon, output_dir,
            split=subont, fold_id=1,
            sparsity_threshold=sparsity_threshold,
        )
        ck.echo(f"[INFO] Wrote {len(written_files)} files for '{subont}'")

        all_summaries.append({
            'subontology':            subont,
            'n_go_terms':             len(all_terms),
            'n_predictions':          len(predictions),
            'n_test_taxa_with_preds': len(predictions_by_taxon),
            'n_output_files':         len(written_files),
        })

        del predictions, predictions_by_taxon, all_terms
        gc.collect()

    # ------------------------------------------------------------------
    # 12. Save summary
    # ------------------------------------------------------------------
    summary = {
        'subontology':           subontology,
        'subontologies_run':     subontologies,
        'n_train_organisms':     len(train_organisms),
        'n_test_organisms':      len(test_organisms),
        'n_train_annotations':   n_train_annotations,
        'evalue':                evalue,
        'max_target_seqs':       max_target_seqs,
        'sparsity_threshold':    sparsity_threshold,
        'output_dir':            output_dir,
        'per_subontology':       all_summaries,
    }
    summary_path = os.path.join(output_dir, 'summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    ck.echo(f"\n[INFO] Summary written to {summary_path}")

    if psutil is not None:
        try:
            mem_gb = psutil.virtual_memory().used / (1024 ** 3)
            ck.echo(f"[INFO] Final memory usage: {mem_gb:.2f} GB")
        except Exception:
            pass

    ck.echo("\n=== DONE ===")
    ck.echo(f"  Sub-ontologies run : {', '.join(subontologies)}")
    ck.echo(f"  Train organisms    : {len(train_organisms)}")
    ck.echo(f"  Test organisms     : {len(test_organisms)}")
    for s in all_summaries:
        ck.echo(f"  [{s['subontology'].upper()}] "
                f"{s['n_predictions']} proteins predicted, "
                f"{s['n_output_files']} files → predictions/{s['subontology']}/")
    ck.echo(f"  Output directory   : {output_dir}")


if __name__ == '__main__':
    main()
