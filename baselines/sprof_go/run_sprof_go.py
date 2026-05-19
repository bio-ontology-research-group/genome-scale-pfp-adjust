"""
Run SPROF-GO predictions for UniProt reference proteomes folds.

Outputs prediction TSV files compatible with pipeline evaluation:
protein_id\tGO:XXXXXX|score\t...
"""

import gzip
import os
import shlex
import subprocess
from typing import Dict, Iterable, List, Optional, Tuple

import click as ck


SUBONTOLOGIES = ['mf', 'bp', 'cc']
SCORE_THRESHOLD = 0.01


def _split_semicolons(line: str) -> List[str]:
    return [tok.strip() for tok in line.split(';') if tok.strip()]


def normalize_protein_id(header: str) -> str:
    """Normalize FASTA header to protein ID."""
    header = header.strip()
    if header.startswith('>'):
        header = header[1:]
    if '|' in header:
        parts = header.split('|')
        if len(parts) >= 3:
            header = parts[2]
    return header.split()[0]


def load_test_organisms(test_organisms_file: str) -> List[str]:
    with open(test_organisms_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def load_proteome_mapping(ids_tsv: str) -> Dict[str, Tuple[str, str]]:
    mapping: Dict[str, Tuple[str, str]] = {}
    with open(ids_tsv, 'r') as f:
        header = next(f, None)
        if header is None:
            return mapping
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            taxon_id, proteome_id, _, _, domain = parts[:5]
            mapping[taxon_id] = (proteome_id, domain)
    return mapping


def find_proteome_fasta_file(main_proteomes_dir: str, proteome_id: str, taxon_id: str) -> Optional[str]:
    domains = ['Bacteria', 'Eukaryota', 'Archaea', 'Viruses']
    for domain in domains:
        proteome_file = os.path.join(
            main_proteomes_dir,
            domain,
            proteome_id,
            f'{proteome_id}_{taxon_id}.fasta.gz',
        )
        if os.path.exists(proteome_file):
            return proteome_file
    return None


def prepare_fasta_for_sprof(proteome_gz: str, output_fasta: str) -> int:
    """Decompress and write FASTA for SPROF-GO."""
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    sequences_written = 0
    with gzip.open(proteome_gz, 'rt') as input_fasta, open(output_fasta, 'w', encoding='utf-8') as output:
        protein_id = ''
        seq_chunks: List[str] = []
        for line in input_fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if protein_id and seq_chunks:
                    output.write(f">{protein_id}\n{''.join(seq_chunks)}\n")
                    sequences_written += 1
                protein_id = normalize_protein_id(line)
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if protein_id and seq_chunks:
            output.write(f">{protein_id}\n{''.join(seq_chunks)}\n")
            sequences_written += 1
    return sequences_written


def run_sprof_go_prediction(
    fasta: str,
    output_dir: str,
    subontology: str,
    sprof_cmd_template: str,
) -> str:
    """Execute SPROF-GO on input FASTA."""
    os.makedirs(output_dir, exist_ok=True)
    formatted_command = sprof_cmd_template.format(
        fasta=fasta,
        output_dir=output_dir.rstrip('/') + '/',
        subontology=subontology,
    )
    command = shlex.split(formatted_command)
    ck.echo(f"[INFO] Running SPROF-GO: {' '.join(command)}")
    subprocess.run(command, check=True)
    return find_sprof_output_file(output_dir)


def split_fasta_file(fasta_path: str, output_dir: str, max_sequences: int) -> List[str]:
    """Split FASTA into chunks with at most max_sequences each."""
    os.makedirs(output_dir, exist_ok=True)
    chunk_paths: List[str] = []
    chunk_index = 0
    sequences_in_chunk = 0
    output_handle = None
    current_chunk_path = None

    def open_new_chunk() -> None:
        nonlocal output_handle, current_chunk_path, chunk_index, sequences_in_chunk
        if output_handle is not None:
            output_handle.close()
        chunk_index += 1
        sequences_in_chunk = 0
        current_chunk_path = os.path.join(output_dir, f"chunk_{chunk_index:04d}.fasta")
        output_handle = open(current_chunk_path, 'w', encoding='utf-8')
        chunk_paths.append(current_chunk_path)

    with open(fasta_path, 'r', encoding='utf-8') as input_fasta:
        for line in input_fasta:
            if line.startswith('>'):
                if output_handle is None or sequences_in_chunk >= max_sequences:
                    open_new_chunk()
                sequences_in_chunk += 1
            if output_handle is None:
                open_new_chunk()
            output_handle.write(line)

    if output_handle is not None:
        output_handle.close()

    return chunk_paths


def find_sprof_output_file(output_dir: str) -> str:
    candidates = [
        os.path.join(output_dir, fname)
        for fname in os.listdir(output_dir)
        if fname.endswith('_all_preds.txt')
    ]
    if not candidates:
        raise FileNotFoundError(f"No *_all_preds.txt files found in {output_dir}")
    candidates.sort(key=os.path.getmtime, reverse=True)
    return candidates[0]


def parse_sprof_go_output(sprof_output: str, subontology: str) -> Dict[str, Dict[str, float]]:
    """Parse SPROF-GO output, return {protein_id: {GO_term: score}}.

    The file begins with a header listing ordered GO term IDs for MF, BP, CC.
    Each protein block contains semicolon-separated scores aligned to those terms.
    """
    if subontology not in SUBONTOLOGIES:
        raise ValueError(f"Unsupported subontology: {subontology}")

    with open(sprof_output, 'r', encoding='utf-8') as fh:
        lines = fh.read().splitlines()
    n = len(lines)
    i = 0

    # --- phase 1: locate and parse header -----------------------------------
    while i < n and 'GO term id and name' not in lines[i]:
        i += 1
    i += 1  # skip the marker line

    global_terms: Dict[str, List[str]] = {}
    for sub in SUBONTOLOGIES:
        label = sub.upper() + ':'
        while i < n and lines[i].strip() != label:
            i += 1
        i += 1  # skip label
        while i < n and not lines[i].strip():
            i += 1
        if i < n:
            global_terms[sub] = [
                t for t in _split_semicolons(lines[i]) if t.startswith('GO:')
            ]
            i += 1
        else:
            global_terms[sub] = []
        # names line (skip)
        while i < n and not lines[i].strip():
            i += 1
        if i < n:
            i += 1

    # --- phase 2: parse per-protein blocks ----------------------------------
    predictions: Dict[str, Dict[str, float]] = {}

    # skip blank lines separating header from first protein
    while i < n and not lines[i].strip():
        i += 1

    while i < n:
        protein_line = lines[i].strip()
        if not protein_line:
            i += 1
            continue
        protein_id = protein_line
        i += 1

        scores_by_sub: Dict[str, List[float]] = {}
        for sub in SUBONTOLOGIES:
            label = sub.upper() + ':'
            while i < n and lines[i].strip() != label:
                if lines[i].strip() and lines[i].strip() not in (
                    'MF:', 'BP:', 'CC:'
                ):
                    # hit next protein before finishing current block
                    break
                i += 1
            if i < n and lines[i].strip() == label:
                i += 1  # skip label
            # skip blank
            while i < n and not lines[i].strip():
                i += 1
            raw_scores: List[float] = []
            if i < n and lines[i].strip() and lines[i].strip() not in (
                'MF:', 'BP:', 'CC:'
            ):
                for tok in _split_semicolons(lines[i]):
                    try:
                        raw_scores.append(float(tok))
                    except ValueError:
                        raw_scores.append(0.0)
                i += 1
            scores_by_sub[sub] = raw_scores

        target_terms = global_terms.get(subontology, [])
        sub_scores = scores_by_sub.get(subontology, [])
        predictions[protein_id] = {
            term: score
            for term, score in zip(target_terms, sub_scores)
            if score >= SCORE_THRESHOLD
        }

        # skip blank lines between proteins
        while i < n and not lines[i].strip():
            i += 1

    return predictions


def merge_predictions(
    predictions_list: Iterable[Dict[str, Dict[str, float]]]
) -> Dict[str, Dict[str, float]]:
    """Merge predictions from multiple runs by taking max score per GO term."""
    merged: Dict[str, Dict[str, float]] = {}
    for predictions in predictions_list:
        for protein_id, terms in predictions.items():
            entry = merged.setdefault(protein_id, {})
            for term, score in terms.items():
                previous = entry.get(term)
                if previous is None or score > previous:
                    entry[term] = score
    return merged


def convert_to_pipeline_format(predictions: Dict[str, Dict[str, float]], output_file: str) -> None:
    """Write predictions in pipeline TSV format."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        for protein_id in sorted(predictions.keys()):
            terms = predictions[protein_id]
            if not terms:
                f.write(f"{protein_id}\n")
                continue
            entries = [f"{term}|{score:.6f}" for term, score in terms.items()]
            f.write(f"{protein_id}\t" + "\t".join(entries) + "\n")


@ck.command()
@ck.option('--test-organisms-file', required=True, type=ck.Path(exists=True, dir_okay=False),
           help='Text file with one test-organism taxon ID per line.')
@ck.option('--proteomes-dir', required=True, type=ck.Path(exists=True, file_okay=False))
@ck.option('--proteomes-ids', required=True, type=ck.Path(exists=True, dir_okay=False))
@ck.option('--output-dir', required=True, type=ck.Path(file_okay=False))
@ck.option('--work-dir', required=True, type=ck.Path(file_okay=False))
@ck.option('--subontology', required=False, default='all', show_default=True, type=ck.Choice(['mf', 'bp', 'cc', 'all'], case_sensitive=False))
@ck.option('--sprof-cmd', required=True, help='Command template for SPROF-GO (use {fasta} and {output_dir} placeholders).')
@ck.option('--fold-id', type=int, default=1, show_default=True,
           help='Fold ID written into predictions_fold_XX_taxon_YYYY.tsv filenames.')
@ck.option('--skip-existing', is_flag=True, help='Skip predictions that already exist.')
@ck.option('--max-taxons', type=int, default=None, help='Optional cap on number of taxons processed.')
@ck.option('--max-seqs', type=int, default=5000, show_default=True, help='Max sequences per SPROF-GO run.')
@ck.option('--dry-run', is_flag=True, help='Print planned actions without executing SPROF-GO.')
def main(
    test_organisms_file: str,
    proteomes_dir: str,
    proteomes_ids: str,
    output_dir: str,
    work_dir: str,
    subontology: str,
    sprof_cmd: str,
    fold_id: int,
    skip_existing: bool,
    max_taxons: Optional[int],
    max_seqs: int,
    dry_run: bool,
) -> None:
    proteome_mapping = load_proteome_mapping(proteomes_ids)
    taxon_ids = load_test_organisms(test_organisms_file)
    if max_taxons:
        taxon_ids = taxon_ids[:max_taxons]
    if not taxon_ids:
        raise ValueError(f"No taxon IDs read from {test_organisms_file}")

    subontologies = ['mf', 'bp', 'cc'] if subontology == 'all' else [subontology]

    for sub in subontologies:
        os.makedirs(os.path.join(output_dir, sub, 'predictions'), exist_ok=True)

    fold_num = fold_id
    ck.echo(f"[INFO] Processing {len(taxon_ids)} test organisms (fold_id={fold_num:02d})")

    for taxon_id in taxon_ids:
        proteome_info = proteome_mapping.get(str(taxon_id))
        if not proteome_info:
            ck.echo(f"[WARNING] No proteome ID found for taxon {taxon_id}")
            continue
        proteome_id, _domain = proteome_info
        proteome_file = find_proteome_fasta_file(proteomes_dir, proteome_id, str(taxon_id))
        if not proteome_file:
            ck.echo(f"[WARNING] Proteome file not found for taxon {taxon_id} (proteome {proteome_id})")
            continue

        fasta_output = os.path.join(work_dir, f"{proteome_id}_{taxon_id}.fasta")
        prediction_outputs = {
            sub: os.path.join(
                output_dir, sub, 'predictions',
                f"predictions_fold_{fold_num:02d}_taxon_{taxon_id}.tsv",
            )
            for sub in subontologies
        }

        if skip_existing and all(os.path.exists(p) for p in prediction_outputs.values()):
            ck.echo(f"[INFO] Skipping existing predictions for taxon {taxon_id}")
            continue

        ck.echo(f"[INFO] Preparing FASTA for taxon {taxon_id} from {proteome_file}")
        if not dry_run:
            sequences_written = prepare_fasta_for_sprof(proteome_file, fasta_output)
            ck.echo(f"[INFO] Wrote {sequences_written} sequences to {fasta_output}")

        sprof_output_dir = os.path.join(work_dir, f"sprof_fold_{fold_num:02d}_taxon_{taxon_id}")
        if dry_run:
            ck.echo(f"[DRY RUN] Would run SPROF-GO for {fasta_output} -> {sprof_output_dir}")
            continue
        if sequences_written > max_seqs:
            ck.echo(f"[INFO] Splitting FASTA into chunks of {max_seqs} sequences")
            chunks_dir = os.path.join(sprof_output_dir, "chunks")
            chunk_fastas = split_fasta_file(fasta_output, chunks_dir, max_seqs)
            print(f"[INFO] Number of chunks: {len(chunk_fastas)}")
            chunk_sprof_outputs = []
            for idx, chunk_fasta in enumerate(chunk_fastas, start=1):
                chunk_output_dir = os.path.join(sprof_output_dir, f"chunk_{idx:04d}")
                sprof_output = run_sprof_go_prediction(
                    chunk_fasta,
                    chunk_output_dir,
                    subontology,
                    sprof_cmd,
                )
                ck.echo(f"[INFO] Collected SPROF-GO output {sprof_output}")
                chunk_sprof_outputs.append(sprof_output)
            for sub in subontologies:
                all_predictions = [parse_sprof_go_output(o, sub) for o in chunk_sprof_outputs]
                predictions = merge_predictions(all_predictions)
                convert_to_pipeline_format(predictions, prediction_outputs[sub])
                ck.echo(f"[INFO] Wrote {sub} predictions to {prediction_outputs[sub]}")
        else:
            sprof_output = run_sprof_go_prediction(fasta_output, sprof_output_dir, subontology, sprof_cmd)
            ck.echo(f"[INFO] Parsing SPROF-GO output {sprof_output}")
            for sub in subontologies:
                predictions = parse_sprof_go_output(sprof_output, sub)
                convert_to_pipeline_format(predictions, prediction_outputs[sub])
                ck.echo(f"[INFO] Wrote {sub} predictions to {prediction_outputs[sub]}")


if __name__ == '__main__':
    main()
