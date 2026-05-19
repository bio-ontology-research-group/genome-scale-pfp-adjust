#!/usr/bin/env python3
"""
Create a train/test/validation split from UniProt proteomes.

Steps:
  1. Sample organisms randomly until cumulative EXPCount >= target_exp_count
     (default 7600, ~10% of 76657). Those organisms form the test set.
  2. From the remaining organisms, collect proteins that carry experimental GO
     annotations and write them to a combined FASTA file.
  3. Cluster all annotated proteins with MMseqs2 easy-cluster.
  4. Split clusters 90/10 (train/validation) so that validation proteins come
     from entirely different clusters than training proteins.

Outputs (all written to --output_dir):
  test_organisms.txt      – one TaxonID per line
  train_organisms.txt     – one TaxonID per line
  annotated_proteins.fasta
  mmseqs_out/             – raw MMseqs2 outputs
  train_proteins.tsv      – protein_id <TAB> taxon_id <TAB> GO:xxx <TAB> ...
  val_proteins.tsv        – same format
  split_info.json         – summary statistics
"""

import argparse
import gzip
import json
import logging
import os
import random
import subprocess
import sys
from collections import defaultdict

import pandas as pd

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Step 1 – Sample test organisms
# ---------------------------------------------------------------------------

def load_consistent_taxons(gaef_file: str) -> set[str]:
    """
    Load taxon IDs marked as Consistent (satisfiable) from a GAEF evaluation file.

    Supports:
    - GAEF eval .out file: Parse "Per-Taxon GAEF Metrics" table; keep taxons
      where the Consistent column (2nd token) is "Yes".
    - per_taxon_gaef_metrics.tsv: Keep rows where satisfiable == "True".

    Returns:
        Set of taxon ID strings.
    """
    if not os.path.exists(gaef_file):
        raise FileNotFoundError(f"GAEF consistent file not found: {gaef_file}")

    consistent: set[str] = set()
    in_table = False

    with open(gaef_file, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if "Per-Taxon GAEF Metrics" in line:
                in_table = True
                continue
            if in_table and "Aggregate GAEF Metrics Summary" in line:
                break
            if in_table and ("---" in line or "Taxon" in line and "Consistent" in line):
                continue
            if in_table and line:
                parts = line.split()
                if len(parts) >= 2:
                    taxon_id = parts[0].strip()
                    consistent_val = parts[1].strip()
                    if consistent_val == "Yes":
                        consistent.add(taxon_id)

    # If no rows found from .out format, try TSV format
    if not consistent:
        with open(gaef_file, "r") as fh:
            lines = fh.readlines()
        if lines:
            header = lines[0].lower().split("\t")
            if "taxon" in header and "satisfiable" in header:
                taxon_idx = header.index("taxon")
                sat_idx = header.index("satisfiable")
                for line in lines[1:]:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) > max(taxon_idx, sat_idx):
                        if parts[sat_idx].strip().lower() == "true":
                            consistent.add(parts[taxon_idx].strip())

    logger.info(
        "Loaded %d consistent taxons from %s",
        len(consistent),
        gaef_file,
    )
    return consistent


def load_proteomes_df(proteomes_ids_file: str) -> pd.DataFrame:
    """Load uniprot_proteomes_ids.tsv into a DataFrame."""
    if not os.path.exists(proteomes_ids_file):
        raise FileNotFoundError(f"Proteomes IDs file not found: {proteomes_ids_file}")
    df = pd.read_csv(proteomes_ids_file, sep="\t")
    required = {"TaxonID", "ProteomeID", "EXPCount", "EntriesCount", "Domain"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in TSV: {missing}")
    df["TaxonID"] = df["TaxonID"].astype(str)
    df["EXPCount"] = pd.to_numeric(df["EXPCount"], errors="coerce").fillna(0).astype(int)
    df["EntriesCount"] = pd.to_numeric(df["EntriesCount"], errors="coerce").fillna(0).astype(int)
    df["Domain"] = df["Domain"].str.lower()
    logger.info(
        "Loaded %d organisms; total EXPCount = %d",
        len(df),
        df["EXPCount"].sum(),
    )
    return df


def sample_test_organisms(df: pd.DataFrame, target_exp_count: int, seed: int):
    """
    Randomly draw organisms (without replacement) until the cumulative EXPCount
    reaches target_exp_count. Returns (test_taxids, remaining_df).
    It tries different random seeds if the initial sample misses the target by more than 100.
    """
    best_taxids = []
    best_cumulative = 0
    best_diff = float('inf')
    
    # Try up to 1000 different shuffles to find a tight bound
    for attempt in range(1000):
        current_seed = seed + attempt
        shuffled = df.sample(frac=1, random_state=current_seed).reset_index(drop=True)
        test_taxids = []
        cumulative = 0
        
        for _, row in shuffled.iterrows():
            if cumulative + row["EXPCount"] > target_exp_count + 100:
                # If adding this pushes us too far over, skip it
                continue
                
            test_taxids.append(row["TaxonID"])
            cumulative += row["EXPCount"]
            
            diff = abs(cumulative - target_exp_count)
            
            # If we're within 100 of the target, we can stop early
            if diff <= 100:
                best_taxids = test_taxids.copy()
                best_cumulative = cumulative
                best_diff = diff
                break
                
        if best_diff <= 100:
            break
            
        # Keep track of the best one we've seen if we don't find a perfect match
        if abs(cumulative - target_exp_count) < best_diff:
            best_diff = abs(cumulative - target_exp_count)
            best_taxids = test_taxids.copy()
            best_cumulative = cumulative

    logger.info(
        "Test set: %d organisms, cumulative EXPCount = %d (target %d)",
        len(best_taxids),
        best_cumulative,
        target_exp_count,
    )
    test_set = set(best_taxids)
    remaining_df = df[~df["TaxonID"].isin(test_set)].copy()
    logger.info("Remaining (train pool): %d organisms", len(remaining_df))
    return best_taxids, remaining_df


# ---------------------------------------------------------------------------
# Step 2 – Collect annotated proteins from remaining organisms
# ---------------------------------------------------------------------------

def find_fasta_gz(main_proteomes_dir: str, proteome_id: str, taxon_id: str) -> str | None:
    """Return the path to a proteome FASTA .gz file, or None if not found."""
    domain_dirs = ["Bacteria", "Eukaryota", "Archaea", "Viruses"]
    for domain in domain_dirs:
        candidate = os.path.join(
            main_proteomes_dir, domain, proteome_id,
            f"{proteome_id}_{taxon_id}.fasta.gz",
        )
        if os.path.exists(candidate):
            return candidate
    return None


def extract_sequences(fasta_gz: str, wanted: set | None = None) -> dict:
    """
    Return {protein_id: sequence} from a .fasta.gz file.

    If *wanted* is given, only sequences whose ID is in that set are kept,
    allowing early termination once all requested proteins have been found.
    """
    sequences = {}
    remaining = set(wanted) if wanted is not None else None
    with gzip.open(fasta_gz, "rt") as fh:
        prot_id = ""
        seq_parts = []
        collect = False
        for line in fh:
            # Stop early when every requested protein has been collected.
            if remaining is not None and not remaining:
                break
            line = line.rstrip()
            if line.startswith(">"):
                if collect and prot_id and seq_parts:
                    sequences[prot_id] = "".join(seq_parts)
                    if remaining is not None:
                        remaining.discard(prot_id)
                # UniProt FASTA header: >sp|ACCESSION|ENTRY_NAME ...
                header = line[1:]
                parts = header.split("|")
                if len(parts) >= 3:
                    prot_id = parts[2].split()[0]
                else:
                    prot_id = header.split()[0]
                collect = (remaining is None) or (prot_id in remaining)
                seq_parts = []
            elif line and collect:
                seq_parts.append(line)
        if collect and prot_id and seq_parts:
            sequences[prot_id] = "".join(seq_parts)
    return sequences


def load_annotations(annotations_dir: str, taxon_id: str) -> dict:
    """
    Load experimental annotations from annots_taxon_{taxon_id}.tsv.
    Returns {protein_id: [go_term, ...]}. Proteins with no annotations
    are excluded.
    """
    annotation_file = os.path.join(annotations_dir, f"annots_taxon_{taxon_id}.tsv")
    if not os.path.exists(annotation_file):
        return {}
    annotations = {}
    with open(annotation_file, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts:
                continue
            protein_id = parts[0].strip()
            go_terms = [t.strip() for t in parts[1:] if t.strip()]
            if protein_id and go_terms:
                annotations[protein_id] = go_terms
    return annotations


def collect_annotated_proteins(
    remaining_df: pd.DataFrame,
    annotations_dir: str,
    main_proteomes_dir: str,
) -> tuple[dict, dict, dict]:
    """
    For every organism in remaining_df, collect proteins that have experimental
    annotations AND have a sequence available.

    Returns:
        protein_annots  : {protein_id: [go_term, ...]}
        protein_seqs    : {protein_id: sequence}
        protein_orgs    : {protein_id: taxon_id}
    """
    protein_annots: dict = {}
    protein_seqs: dict = {}
    protein_orgs: dict = {}

    total_orgs = len(remaining_df)
    for idx, row in enumerate(remaining_df.itertuples(), 1):
        taxon_id = str(row.TaxonID)
        proteome_id = str(row.ProteomeID)

        if idx % 50 == 0 or idx == total_orgs:
            logger.info("  Processing organism %d/%d (taxon %s)", idx, total_orgs, taxon_id)

        annots = load_annotations(annotations_dir, taxon_id)
        if not annots:
            continue

        fasta_gz = find_fasta_gz(main_proteomes_dir, proteome_id, taxon_id)
        if fasta_gz is None:
            logger.debug("No FASTA found for taxon %s / proteome %s", taxon_id, proteome_id)
            continue

        seqs = extract_sequences(fasta_gz, wanted=set(annots.keys()))

        for prot_id, go_terms in annots.items():
            if prot_id in seqs:
                protein_annots[prot_id] = go_terms
                protein_seqs[prot_id] = seqs[prot_id]
                protein_orgs[prot_id] = taxon_id

    logger.info(
        "Collected %d annotated proteins with sequences from %d organisms",
        len(protein_annots),
        len(set(protein_orgs.values())),
    )
    return protein_annots, protein_seqs, protein_orgs


def write_fasta(protein_seqs: dict, output_fasta: str) -> None:
    """Write sequences to a FASTA file."""
    with open(output_fasta, "w") as fh:
        for prot_id, seq in protein_seqs.items():
            fh.write(f">{prot_id}\n{seq}\n")
    logger.info("Wrote %d sequences to %s", len(protein_seqs), output_fasta)


# ---------------------------------------------------------------------------
# Step 3 – MMseqs2 clustering
# ---------------------------------------------------------------------------

def run_mmseqs_cluster(
    fasta_file: str,
    output_dir: str,
    min_seq_id: float = 0.3,
    coverage: float = 0.8,
    threads: int = 1,
) -> str:
    """
    Run `mmseqs easy-cluster` and return the path to the cluster TSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, "mmseqs_clusters")
    tmp_dir = os.path.join(output_dir, "mmseqs_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    cmd = [
        "mmseqs", "easy-cluster",
        fasta_file,
        prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", "0",
        "--threads", str(threads),
        "-v", "1",
    ]
    logger.info("Running MMseqs2: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error("MMseqs2 stdout:\n%s", result.stdout)
        logger.error("MMseqs2 stderr:\n%s", result.stderr)
        raise RuntimeError(f"MMseqs2 easy-cluster failed (exit {result.returncode})")

    cluster_tsv = f"{prefix}_cluster.tsv"
    if not os.path.exists(cluster_tsv):
        raise FileNotFoundError(f"Expected cluster TSV not found: {cluster_tsv}")
    logger.info("MMseqs2 clustering complete. Cluster TSV: %s", cluster_tsv)
    return cluster_tsv


def parse_cluster_tsv(cluster_tsv: str) -> dict:
    """
    Parse MMseqs2 cluster TSV (rep_id <TAB> member_id) into
    {rep_id: [member_id, ...]} where the rep is also included in its cluster.
    """
    clusters: dict = defaultdict(list)
    with open(cluster_tsv, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters[rep].append(member)
    logger.info("Parsed %d clusters from %s", len(clusters), cluster_tsv)
    return dict(clusters)


# ---------------------------------------------------------------------------
# Step 4 – Cluster-based 90/10 train/validation split
# ---------------------------------------------------------------------------

def split_clusters(
    clusters: dict,
    total_proteins: int,
    train_ratio: float = 0.9,
    seed: int = 42,
) -> tuple[set, set]:
    """
    Assign clusters to train or validation sets using greedy bin-packing so
    that the validation set contains approximately (1 - train_ratio) of all
    proteins and no protein appears in both splits.

    Returns (train_protein_ids, val_protein_ids).
    """
    val_target = total_proteins * (1.0 - train_ratio)

    # Sort clusters largest-first for a tighter greedy approximation, then
    # shuffle ties randomly to avoid systematic bias toward any organism.
    rng = random.Random(seed)
    cluster_items = list(clusters.items())
    rng.shuffle(cluster_items)
    cluster_items.sort(key=lambda x: len(x[1]), reverse=True)

    train_proteins: set = set()
    val_proteins: set = set()
    val_count = 0

    for rep, members in cluster_items:
        member_set = set(members)
        if val_count < val_target:
            val_proteins.update(member_set)
            val_count += len(member_set)
        else:
            train_proteins.update(member_set)

    logger.info(
        "Split result: %d train proteins / %d val proteins (target val=%.0f)",
        len(train_proteins),
        len(val_proteins),
        val_target,
    )
    return train_proteins, val_proteins


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_protein_tsv(
    protein_ids: set,
    protein_annots: dict,
    protein_orgs: dict,
    output_file: str,
) -> int:
    """
    Write a TSV: protein_id <TAB> taxon_id <TAB> GO:xxx <TAB> GO:yyy ...
    Returns the number of proteins written.
    """
    written = 0
    with open(output_file, "w") as fh:
        for prot_id in sorted(protein_ids):
            taxon_id = protein_orgs.get(prot_id, "")
            go_terms = protein_annots.get(prot_id, [])
            line = "\t".join([prot_id, taxon_id] + go_terms)
            fh.write(line + "\n")
            written += 1
    logger.info("Wrote %d proteins to %s", written, output_file)
    return written


def write_organism_list(taxon_ids: list, output_file: str) -> None:
    with open(output_file, "w") as fh:
        for tid in taxon_ids:
            fh.write(f"{tid}\n")
    logger.info("Wrote %d taxon IDs to %s", len(taxon_ids), output_file)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create train/test/validation split from UniProt proteomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--proteomes_ids_file",
        default="data/uniprot_proteomes_ids.tsv",
        help="Path to uniprot_proteomes_ids.tsv",
    )
    parser.add_argument(
        "--annotations_dir",
        default="${DATA_DIR}/swissprot_proteomes_folds/annotations-go-basic",
        help="Directory containing annots_taxon_<taxid>.tsv files",
    )
    parser.add_argument(
        "--main_proteomes_dir",
        default="${DATA_DIR}/uniprot_reference_proteomes",
        help="Base directory with UniProt reference proteome FASTA files",
    )
    parser.add_argument(
        "--output_dir",
        default="splits/heldout",
        help="Directory to write all outputs",
    )
    parser.add_argument(
        "--target_exp_count",
        type=int,
        default=7600,
        help="Cumulative EXPCount target for the test set (~10%% of total)",
    )
    parser.add_argument(
        "--min_seq_id",
        type=float,
        default=0.3,
        help="MMseqs2 minimum sequence identity for clustering",
    )
    parser.add_argument(
        "--train_ratio",
        type=float,
        default=0.9,
        help="Fraction of annotated proteins assigned to training",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=os.cpu_count() or 1,
        help="Number of threads for MMseqs2",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--gaef_consistent_file",
        default=None,
        help="Path to GAEF eval .out or per_taxon_gaef_metrics.tsv. If set, only taxons marked Consistent (satisfiable) are eligible for the test set.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Add a file handler now that we know the output dir
    file_handler = logging.FileHandler(
        os.path.join(args.output_dir, "create_heldout_split.log"), mode="w"
    )
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Configuration ===")
    for k, v in vars(args).items():
        logger.info("  %s = %s", k, v)

    # ------------------------------------------------------------------
    # Step 1: Sample test organisms
    # ------------------------------------------------------------------
    logger.info("\n=== Step 1: Sampling test organisms ===")
    df = load_proteomes_df(args.proteomes_ids_file)

    gaef_filter_info = None
    if args.gaef_consistent_file:
        consistent_taxons = load_consistent_taxons(args.gaef_consistent_file)
        df_test_pool = df[df["TaxonID"].isin(consistent_taxons)].copy()
        test_pool_exp = int(df_test_pool["EXPCount"].sum())
        gaef_filter_info = {
            "file": args.gaef_consistent_file,
            "n_consistent_taxons": len(consistent_taxons),
            "test_pool_exp_count": test_pool_exp,
        }
        logger.info(
            "GAEF filter: %d consistent organisms in test pool (EXPCount %d)",
            len(df_test_pool),
            test_pool_exp,
        )
    else:
        df_test_pool = df

    test_taxids, _ = sample_test_organisms(
        df_test_pool, args.target_exp_count, args.seed
    )
    remaining_df = df[~df["TaxonID"].isin(test_taxids)].copy()
    train_taxids = remaining_df["TaxonID"].tolist()

    test_org_file = os.path.join(args.output_dir, "test_organisms.txt")
    train_org_file = os.path.join(args.output_dir, "train_organisms.txt")
    write_organism_list(test_taxids, test_org_file)
    write_organism_list(train_taxids, train_org_file)

    # ------------------------------------------------------------------
    # Step 2: Collect annotated proteins from remaining organisms
    # ------------------------------------------------------------------
    logger.info("\n=== Step 2: Collecting annotated proteins ===")
    protein_annots, protein_seqs, protein_orgs = collect_annotated_proteins(
        remaining_df, args.annotations_dir, args.main_proteomes_dir
    )

    if not protein_annots:
        logger.error(
            "No annotated proteins found. Check --annotations_dir and "
            "--main_proteomes_dir paths."
        )
        sys.exit(1)

    fasta_path = os.path.join(args.output_dir, "annotated_proteins.fasta")
    write_fasta(protein_seqs, fasta_path)

    # ------------------------------------------------------------------
    # Step 3: MMseqs2 clustering
    # ------------------------------------------------------------------
    logger.info("\n=== Step 3: Clustering with MMseqs2 ===")
    mmseqs_dir = os.path.join(args.output_dir, "mmseqs_out")
    cluster_tsv = run_mmseqs_cluster(
        fasta_path,
        mmseqs_dir,
        min_seq_id=args.min_seq_id,
        coverage=0.8,
        threads=args.threads,
    )
    clusters = parse_cluster_tsv(cluster_tsv)
    print(f"\nNumber of clusters found by MMseqs2: {len(clusters)}\n")

    # ------------------------------------------------------------------
    # Step 4: Cluster-based 90/10 train/val split
    # ------------------------------------------------------------------
    logger.info("\n=== Step 4: Splitting train/validation by cluster ===")
    total_proteins = len(protein_annots)
    train_proteins, val_proteins = split_clusters(
        clusters,
        total_proteins,
        train_ratio=args.train_ratio,
        seed=args.seed,
    )

    # Restrict to proteins that have annotations (in case clustering added
    # proteins not in our annotated set – shouldn't happen but be safe)
    train_proteins &= set(protein_annots.keys())
    val_proteins &= set(protein_annots.keys())

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    logger.info("\n=== Writing outputs ===")
    train_tsv = os.path.join(args.output_dir, "train_proteins.tsv")
    val_tsv = os.path.join(args.output_dir, "val_proteins.tsv")
    n_train = write_protein_tsv(train_proteins, protein_annots, protein_orgs, train_tsv)
    n_val = write_protein_tsv(val_proteins, protein_annots, protein_orgs, val_tsv)

    # Summary
    test_exp_count = df[df["TaxonID"].isin(test_taxids)]["EXPCount"].sum()
    train_exp_count = remaining_df["EXPCount"].sum()

    split_info = {
        "seed": args.seed,
        "target_exp_count": args.target_exp_count,
        "min_seq_id": args.min_seq_id,
        "train_ratio": args.train_ratio,
        "test_set": {
            "n_organisms": len(test_taxids),
            "cumulative_exp_count": int(test_exp_count),
            "taxon_ids": test_taxids,
        },
        "train_pool": {
            "n_organisms": len(train_taxids),
            "cumulative_exp_count": int(train_exp_count),
        },
        "annotated_proteins": {
            "total": total_proteins,
            "n_clusters": len(clusters),
        },
        "train_split": {
            "n_proteins": n_train,
        },
        "val_split": {
            "n_proteins": n_val,
        },
    }
    if gaef_filter_info is not None:
        split_info["gaef_consistent_filter"] = gaef_filter_info

    split_info_path = os.path.join(args.output_dir, "split_info.json")
    with open(split_info_path, "w") as fh:
        json.dump(split_info, fh, indent=2)
    logger.info("Wrote split summary to %s", split_info_path)

    logger.info("\n=== Done ===")
    logger.info("  Test organisms:          %d  (EXPCount %d)", len(test_taxids), test_exp_count)
    logger.info("  Train-pool organisms:    %d  (EXPCount %d)", len(train_taxids), train_exp_count)
    logger.info("  Annotated proteins:      %d", total_proteins)
    logger.info("  Clusters:                %d", len(clusters))
    logger.info("  Train proteins:          %d  (%.1f%%)", n_train, 100 * n_train / total_proteins)
    logger.info("  Validation proteins:     %d  (%.1f%%)", n_val, 100 * n_val / total_proteins)


if __name__ == "__main__":
    main()
