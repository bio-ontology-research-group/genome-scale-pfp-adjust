#!/usr/bin/env python3
"""
Regenerate the per-protein train/val TSVs for the heldout split deterministically.

create_heldout_split.py samples the test organisms and then builds
train_proteins.tsv / val_proteins.tsv. Only the organism lists and split_info.json
were preserved with the manuscript, not the per-protein TSVs. This wrapper
reproduces those TSVs from the *preserved* organism lists, so the test set is held
fixed and only the train/val partition is rebuilt. With MMseqs2 run single-threaded
and the cluster split seeded (random.Random(42)), the result is deterministic.

Steps (Steps 2-4 of create_heldout_split.py; Step 1 is skipped because the test
set is read from the committed test_organisms.txt):
  1. read train_organisms.txt -> remaining_df
  2. collect annotated proteins with sequences
  3. MMseqs2 easy-cluster (threads=1)
  4. cluster-level 90/10 train/val split (seed 42)

Usage (on a node with `mmseqs` on PATH):
    python data_prep/regen_heldout_proteins.py \
        --split_dir splits/heldout \
        --proteomes_ids_file data/uniprot_proteomes_ids.tsv \
        --annotations_dir ${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic \
        --main_proteomes_dir ${DATA_DIR}/uniprot_reference_proteomes \
        --threads 1 --seed 42
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from data_prep.create_heldout_split import (  # noqa: E402
    load_proteomes_df,
    collect_annotated_proteins,
    write_fasta,
    run_mmseqs_cluster,
    parse_cluster_tsv,
    split_clusters,
    write_protein_tsv,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--split_dir", default="splits/heldout",
                    help="Dir with preserved train_organisms.txt / test_organisms.txt; "
                         "train_proteins.tsv / val_proteins.tsv are written here.")
    ap.add_argument("--proteomes_ids_file", default="data/uniprot_proteomes_ids.tsv")
    ap.add_argument("--annotations_dir", required=True)
    ap.add_argument("--main_proteomes_dir", required=True)
    ap.add_argument("--train_ratio", type=float, default=0.9)
    ap.add_argument("--threads", type=int, default=1,
                    help="MMseqs2 threads; 1 for bit-reproducible clustering.")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    train_org_file = os.path.join(args.split_dir, "train_organisms.txt")
    with open(train_org_file) as fh:
        train_taxids = {ln.strip() for ln in fh if ln.strip()}
    print(f"[INFO] {len(train_taxids)} preserved train organisms from {train_org_file}")

    df = load_proteomes_df(args.proteomes_ids_file)
    remaining_df = df[df["TaxonID"].isin(train_taxids)].copy()
    print(f"[INFO] matched {len(remaining_df)} organisms in proteomes table")

    protein_annots, protein_seqs, protein_orgs = collect_annotated_proteins(
        remaining_df, args.annotations_dir, args.main_proteomes_dir)
    if not protein_annots:
        sys.exit("[ERROR] no annotated proteins collected; check data paths")

    fasta_path = os.path.join(args.split_dir, "annotated_proteins.fasta")
    write_fasta(protein_seqs, fasta_path)

    mmseqs_dir = os.path.join(args.split_dir, "mmseqs_out")
    cluster_tsv = run_mmseqs_cluster(fasta_path, mmseqs_dir, min_seq_id=0.3,
                                     coverage=0.8, threads=args.threads)
    clusters = parse_cluster_tsv(cluster_tsv)

    total_proteins = len(protein_annots)
    train_proteins, val_proteins = split_clusters(
        clusters, total_proteins, train_ratio=args.train_ratio, seed=args.seed)
    train_proteins &= set(protein_annots.keys())
    val_proteins &= set(protein_annots.keys())

    train_tsv = os.path.join(args.split_dir, "train_proteins.tsv")
    val_tsv = os.path.join(args.split_dir, "val_proteins.tsv")
    n_train = write_protein_tsv(train_proteins, protein_annots, protein_orgs, train_tsv)
    n_val = write_protein_tsv(val_proteins, protein_annots, protein_orgs, val_tsv)

    # Cross-check against the committed split_info.json counts (informational).
    info_path = os.path.join(args.split_dir, "split_info.json")
    if os.path.exists(info_path):
        with open(info_path) as fh:
            info = json.load(fh)
        print(f"[INFO] regenerated train={n_train} val={n_val}; "
              f"committed train={info.get('train_split', {}).get('n_proteins')} "
              f"val={info.get('val_split', {}).get('n_proteins')} "
              f"clusters now={len(clusters)} committed={info.get('annotated_proteins', {}).get('n_clusters')}")
    print(f"[INFO] wrote {train_tsv} and {val_tsv}")


if __name__ == "__main__":
    main()
