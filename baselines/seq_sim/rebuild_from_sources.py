"""Regenerate sequence-transfer scores from verified source labels and sequences."""
import argparse
from contextlib import ExitStack
import gzip
import hashlib
from itertools import groupby
import json
import math
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from data_prep.prepare_baseline_inputs import digest
from gaef_patches.ontology import Ontology, NAMESPACES


def sequences(path):
    protein = None
    sequence = []
    in_sequence = False
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('ID   '):
                if protein is not None:
                    raise ValueError('Unterminated DAT entry')
                protein = line.split()[1]
                sequence = []
            elif line.startswith('SQ   '):
                in_sequence = True
            elif line.startswith('//'):
                if not protein or not sequence:
                    raise ValueError('Missing protein ID/sequence')
                yield protein, ''.join(sequence)
                protein = None
                in_sequence = False
            elif in_sequence:
                text = ''.join(line.split())
                if not text.isalpha():
                    raise ValueError(f'Invalid sequence for {protein}')
                sequence.append(text)
    if protein is not None:
        raise ValueError('Truncated DAT file')


def transfer_scores(hits, annotations):
    total = 0.0
    sums = {}
    for subject, bitscore in hits.items():
        if subject not in annotations or not math.isfinite(bitscore) or bitscore <= 0:
            raise ValueError(f'Invalid alignment subject/score: {subject}')
        total += bitscore
        for term in annotations[subject]:
            sums[term] = sums.get(term, 0.0) + bitscore
    return {term:value/total for term,value in sums.items()} if total else {}


def run(args):
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Use a new output directory')
    refs = Path(args.references)
    manifest_path = refs/'reference-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest['go_sha256'] != digest(args.go):
        raise ValueError('Reference ontology differs')
    records = {r['taxon']:r for r in manifest['records']}
    train = set(Path(args.train_taxa).read_text().split())
    test = set(Path(args.test_taxa).read_text().split())
    if not train or not test or train & test or not test <= records.keys():
        raise ValueError('Empty, overlapping or missing organism sets')
    go = Ontology(args.go)
    terms = {}
    inputs = {str(p):digest(p) for p in [manifest_path,Path(args.go),Path(args.train_taxa),Path(args.test_taxa)]}
    for aspect in NAMESPACES:
        path = Path(args.terms_root)/aspect/'terms.pkl'
        inputs[str(path)] = digest(path)
        terms[aspect] = set(pd.read_pickle(path)['gos'])
        if any(go.canonical_id(t) != t or go.get_namespace(t) != NAMESPACES[aspect] for t in terms[aspect]):
            raise ValueError('Vocabulary is not canonical or has another aspect')
    labels = {}
    test_proteins = {}
    by_taxon = {}
    for taxon in sorted((train | test) & records.keys()):
        path = refs/records[taxon]['reference_file']
        if digest(path) != records[taxon]['reference_sha256']:
            raise ValueError(f'Reference checksum differs: {path}')
        ids = set()
        with path.open() as handle:
            for line in handle:
                protein,*annotations = line.rstrip('\n').split('\t')
                if not protein or protein in ids:
                    raise ValueError('Empty or duplicate reference protein')
                ids.add(protein)
                if taxon in test:
                    if protein in test_proteins:
                        raise ValueError('Protein belongs to multiple test taxa')
                    test_proteins[protein] = taxon
                else:
                    annotations = {t for t in annotations if t}
                    if annotations:
                        if protein in labels:
                            raise ValueError('Protein belongs to multiple training taxa')
                        labels[protein] = annotations
        by_taxon[taxon] = ids
    if set(labels) & test_proteins.keys():
        raise ValueError('Training/query protein overlap')
    out.mkdir(parents=True,exist_ok=True)
    tmp = out/'scratch';tmp.mkdir()
    train_fasta, query_fasta = out/'train.fasta',out/'query.fasta'
    found_train = set()
    with train_fasta.open('w') as train_handle, query_fasta.open('w') as query_handle:
        for taxon, expected in by_taxon.items():
            row = records[taxon]
            path = Path(args.source_root).joinpath(*Path(row['source_file']).parts[-3:])
            if digest(path) != row['source_sha256']:
                raise ValueError(f'Source checksum differs: {path}')
            inputs[str(path)] = row['source_sha256']
            seen = set()
            for protein, sequence in sequences(path):
                if protein in seen:
                    raise ValueError('Duplicate source entry')
                seen.add(protein)
                if taxon in test:
                    query_handle.write(f'>{protein}\n{sequence}\n')
                elif protein in labels:
                    train_handle.write(f'>{protein}\n{sequence}\n')
                    found_train.add(protein)
            if seen != expected:
                raise ValueError(f'Source/reference protein universe differs for {taxon}')
    if found_train != labels.keys():
        raise ValueError('Training labels have no sequence')
    version = subprocess.check_output([args.diamond,'version'],text=True).strip()
    database = out/'train'
    commands = [[args.diamond,'makedb','--in',str(train_fasta),'-d',str(database),'--threads',str(args.threads)],
        [args.diamond,'blastp','-q',str(query_fasta),'-d',str(database),'-o',str(out/'alignments.tsv'),
         '--outfmt','6','qseqid','sseqid','bitscore','evalue','--evalue','1e-5','--max-target-seqs','50',
         '--threads',str(args.threads),'--sensitive','--block-size','2.0','--index-chunks','4',
         '--min-score','50','--gapopen','11','--gapextend','1','--matrix','BLOSUM62','-t',str(tmp)]]
    for command in commands:
        subprocess.run(command,check=True)
    output_files = {}
    with ExitStack() as stack:
        handles = {}
        for aspect in terms:
            folder = out/aspect/'predictions';folder.mkdir(parents=True)
            for taxon in sorted(test):
                path = folder/f'predictions_fold_01_taxon_{taxon}.tsv'
                handles[aspect,taxon] = stack.enter_context(path.open('w'))
                output_files[aspect,taxon] = path
        seen = set()
        def write(protein,scores):
            taxon = test_proteins[protein]
            for aspect, vocabulary in terms.items():
                tokens = [f'{term}|{value:.17g}' for term,value in sorted(scores.items())
                          if term in vocabulary and value >= args.score_floor]
                handles[aspect,taxon].write('\t'.join([protein]+tokens)+'\n')
        with (out/'alignments.tsv').open() as alignments:
            rows = (line.rstrip('\n').split('\t') for line in alignments)
            for protein, group in groupby(rows,key=lambda row:row[0]):
                if protein not in test_proteins or protein in seen:
                    raise ValueError('Unknown or non-contiguous alignment query')
                hits = {}
                for row in group:
                    if len(row) != 4:
                        raise ValueError('Malformed alignment')
                    subject, score = row[1],float(row[2])
                    hits[subject] = max(hits.get(subject,0),score)
                write(protein,transfer_scores(hits,labels));seen.add(protein)
        for protein in sorted(test_proteins.keys()-seen):
            write(protein,{})
    result = dict(schema_version=1,inputs=inputs,reference_manifest_sha256=digest(manifest_path),
                  train_taxa=len(train),available_train_taxa=len(train & records.keys()),
                  unavailable_train_taxa=sorted(train-records.keys()),test_taxa=len(test),
                  training_proteins=len(labels),query_proteins=len(test_proteins),queries_with_hits=len(seen),
                  score_floor=args.score_floor,diamond_version=version,commands=commands,
                  outputs={str(p.relative_to(out)):digest(p) for p in output_files.values()},
                  derived_inputs={p.name:digest(p) for p in [train_fasta,query_fasta,out/'alignments.tsv',out/'train.dmnd']},
                  code_sha256=digest(__file__))
    (out/'prediction-manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs','outputs','derived_inputs']},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['references','source-root','go','terms-root','train-taxa','test-taxa','output']:
        parser.add_argument('--'+name,required=True)
    parser.add_argument('--diamond',default='diamond')
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--score-floor',type=float,default=.001)
    args=parser.parse_args()
    if args.threads < 1 or not 0 <= args.score_floor <= .001:
        parser.error('Positive thread count and score floor in [0, 0.001] required')
    run(args)
