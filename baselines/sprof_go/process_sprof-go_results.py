"""
Process SPROF-GO results to convert them to TSV format.

Under ${DATA_DIR}/sprof_go/sprof_go_work/sprof_fold_XX_taxon_YYYY/chunk_XXXX/chunk_XXXX_all_preds.txt, we have the SPROF-GO predictions written in this format:
```
GO term id and name:
MF:
GO:0000030; GO:0000049; ...
<name1>; <name2>; ...
BP:
<GO:XXXXX>; <GO:YYYYY>; ...
<name1>; <name2>; ...
CC:
<GO:XXXXX>; <GO:YYYYY>; ...
<name1>; <name2>; ...

<protein_id1>
MF:
0.0; 0.001; ...
BP:
0.0; 0.001; ...
CC:
0.006; 0.007; ...

<protein_id2>
...
```

This script converts the above format to the following TSV format:
```
protein_id\tGO:XXXXX|score\tGO:YYYYY|score...
protein_id\tGO:XXXXX|score\tGO:YYYYY|score...
...
```
with one file per taxon and subontology combining all chunks.

Output files into ${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results/<subontology>/predictions/

"""

import concurrent.futures
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

import click as ck
import time

sys.path.insert(0, os.path.dirname(__file__))
from run_sprof_go import merge_predictions, convert_to_pipeline_format

SUBONTOLOGIES = ['mf', 'bp', 'cc']
_FOLD_TAXON_RE = re.compile(r'sprof_fold_(\d+)_taxon_(\w+)')


def _split_semicolons(line: str) -> List[str]:
    return [tok.strip() for tok in line.split(';') if tok.strip()]


def parse_sprof_go_file_all(
    filepath: str,
    subs: List[str],
    score_threshold: float,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Parse a SPROF-GO *_all_preds.txt file for all requested subontologies in one pass.

    The file begins with a header that lists ordered GO term IDs for MF, BP and
    CC. Each subsequent protein block contains semicolon-separated scores that
    align positionally to those global term lists.

    Returns {sub: {protein_id: {go_term: score}}}.
    """
    with open(filepath, 'r', encoding='utf-8') as fh:
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
    predictions: Dict[str, Dict[str, Dict[str, float]]] = {sub: {} for sub in subs}

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

        for sub in subs:
            target_terms = global_terms.get(sub, [])
            sub_scores = scores_by_sub.get(sub, [])
            predictions[sub][protein_id] = {
                term: score
                for term, score in zip(target_terms, sub_scores)
                if score >= score_threshold
            }

        # skip blank lines between proteins
        while i < n and not lines[i].strip():
            i += 1

    return predictions


def parse_sprof_go_file(
    filepath: str,
    subontology: str,
    score_threshold: float,
) -> Dict[str, Dict[str, float]]:
    """Parse a SPROF-GO *_all_preds.txt file for one subontology.

    Thin wrapper around parse_sprof_go_file_all for single-subontology use.
    """
    return parse_sprof_go_file_all(filepath, [subontology], score_threshold)[subontology]


def _parse_chunk_worker(
    args: Tuple[str, List[str], float],
) -> Dict[str, Dict[str, Dict[str, float]]]:
    filepath, subs, score_threshold = args
    return parse_sprof_go_file_all(filepath, subs, score_threshold)


def collect_chunk_files(work_dir: str) -> Dict[Tuple[str, str], List[str]]:
    """Walk work_dir and group *_all_preds.txt paths by (fold_id, taxon_id)."""
    groups: Dict[Tuple[str, str], List[str]] = defaultdict(list)
    for dirpath, _dirnames, filenames in os.walk(work_dir):
        for fname in filenames:
            if not fname.endswith('_all_preds.txt'):
                continue
            full_path = os.path.join(dirpath, fname)
            m = _FOLD_TAXON_RE.search(dirpath)
            if not m:
                ck.echo(f"[WARNING] Could not parse fold/taxon from path: {dirpath}", err=True)
                continue
            fold_id, taxon_id = m.group(1), m.group(2)
            groups[(fold_id, taxon_id)].append(full_path)
    return groups


@ck.command()
@ck.option('--work-dir', default='${DATA_DIR}/sprof_go/sprof_go_work',
           type=ck.Path(exists=True, file_okay=False), show_default=True,
           help='Root of the sprof_go_work directory containing sprof_fold_XX_taxon_YY subdirs.')
@ck.option('--output-dir', default='${DATA_DIR}/swissprot_proteomes_folds/sprof-go_results',
           type=ck.Path(file_okay=False), show_default=True,
           help='Root output directory; results go to <output-dir>/<subontology>/predictions/.')
@ck.option('--subontology', default='all', show_default=True,
           type=ck.Choice(['mf', 'bp', 'cc', 'all'], case_sensitive=False),
           help='Subontology to process.')
@ck.option('--score-threshold', default=0.01, show_default=True, type=float,
           help='Minimum score to include a GO term.')
@ck.option('--skip-existing', is_flag=True,
           help='Skip writing output files that already exist.')
@ck.option('--workers', default=0, show_default=True, type=int,
           help='Number of parallel worker processes per taxon group. 0 = os.cpu_count().')
def main(
    work_dir: str,
    output_dir: str,
    subontology: str,
    score_threshold: float,
    skip_existing: bool,
    workers: int,
) -> None:
    subs = SUBONTOLOGIES if subontology == 'all' else [subontology]
    start_time = time.time()
    for sub in subs:
        os.makedirs(os.path.join(output_dir, sub, 'predictions'), exist_ok=True)

    groups = collect_chunk_files(work_dir)
    if not groups:
        ck.echo('[ERROR] No *_all_preds.txt files found under work_dir.', err=True)
        raise SystemExit(1)

    ck.echo(f"[INFO] Found {len(groups)} (fold, taxon) groups across {sum(len(v) for v in groups.values())} chunk files.")

    for (fold_id, taxon_id), chunk_files in sorted(groups.items()):
        fold_num = int(fold_id)
        ck.echo(f"[INFO] Processing fold {fold_num:02d} taxon {taxon_id} ({len(chunk_files)} chunk(s))")
        chunk_start_time = time.time()

        out_paths = {
            sub: os.path.join(
                output_dir, sub, 'predictions',
                f"predictions_fold_{fold_num:02d}_taxon_{taxon_id}.tsv",
            )
            for sub in subs
        }

        # Determine which subs still need writing before touching the files
        subs_needed = []
        for sub in subs:
            if skip_existing and os.path.exists(out_paths[sub]):
                ck.echo(f"[INFO] Skipping existing {out_paths[sub]}")
            else:
                subs_needed.append(sub)

        if not subs_needed:
            continue

        # Parse each chunk file once for all needed subontologies
        n_workers = workers or os.cpu_count()
        args_list = [(cf, subs_needed, score_threshold) for cf in sorted(chunk_files)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(n_workers, len(args_list))) as pool:
            chunk_results = list(pool.map(_parse_chunk_worker, args_list))

        for sub in subs_needed:
            all_preds = [chunk[sub] for chunk in chunk_results]
            merged = merge_predictions(all_preds)
            convert_to_pipeline_format(merged, out_paths[sub])
            ck.echo(f"[INFO] Wrote {len(merged)} proteins ({sub}) -> {out_paths[sub]}")
        chunk_end_time = time.time()
        ck.echo(f"[INFO] Time taken for fold {fold_num:02d} taxon {taxon_id}: {chunk_end_time - chunk_start_time:.2f} seconds")
    end_time = time.time()
    ck.echo(f"[INFO] Total time taken: {end_time - start_time:.2f} seconds")


if __name__ == '__main__':
    main()