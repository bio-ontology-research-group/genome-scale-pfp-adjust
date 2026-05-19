#!/usr/bin/env python3
"""
Build a reproducible dataset:
1) Download Swiss-Prot flatfile (uniprot_sprot.dat.gz)
2) Extract proteins with experimental GO annotations (EXP|IDA|IPI|IMP|IGI|IEP|TAS|IC|HTP|HDA|HMP|HGI|HEP)
3) Collect their organisms (NCBI TaxIDs) from OX lines
4) Download genome protein FASTA data for those organisms from NCBI using the `datasets` CLI
5) Write checksums + metadata for provenance

Usage:
  python build_dataset.py --out ./my_dataset
Optional:
  --max-taxa 100          # only download NCBI proteomes for the first N taxa
  --taxid-include FILE    # newline-delimited taxids to *limit to* (intersect with Swiss-Prot-derived set)
  --taxid-exclude FILE    # newline-delimited taxids to exclude
  --refseq-only           # prefer RefSeq assemblies when using NCBI Datasets
  --organism-report       # also emit organism names in the taxon list
  --threads 4             # parallel downloads from NCBI (default 2)
  --skip-ncbi             # only prepare Swiss-Prot-derived taxon list, skip NCBI downloads
  --swissprot-url URL     # override Swiss-Prot URL if needed

Requirements:
  - Python 3.9+
  - tqdm (pip install tqdm)
  - NCBI Datasets CLI in PATH (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/), command name: `datasets`

This script is idempotent; existing downloads are skipped when checksums match.
"""
import argparse
import concurrent.futures as cf
import datetime
import gzip
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
import zipfile
from typing import Dict, List, Tuple
from urllib.request import urlopen, Request

SWISSPROT_DEFAULT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz"

def sha256_of_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def download(url: str, dst: str, chunk_size: int = 1 << 20):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + ".part"
    req = Request(url, headers={"User-Agent": "dataset-builder/1.0"})
    with urlopen(req) as r, open(tmp, "wb") as f:
        while True:
            b = r.read(chunk_size)
            if not b: break
            f.write(b)
    os.replace(tmp, dst)

def parse_uniprot_dat_for_taxa(dat_gz_path: str, want_names: bool = False) -> Tuple[Dict[str,int], Dict[str,str]]:
    """
    Parse Swiss-Prot .dat.gz, collecting entries with experimental GO annotations and extracting NCBI TaxID from OX line.
    Returns:
      counts: taxid -> number of proteins with experimental GO annotations
      names:  taxid -> organism name (from OS)  [only if want_names]
    """
    tax_counts: Dict[str,int] = {}
    tax_names: Dict[str,str] = {}

    curr_has_exp_go = False
    curr_taxid = None
    curr_os = None

    # Experimental evidence codes for GO annotations
    EXP_CODES = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}

    ox_tax_re = re.compile(r"NCBI_TaxID=(\d+)")
    with gzip.open(dat_gz_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID "):
                # start of entry
                curr_has_exp_go = False
                curr_taxid = None
                curr_os = None
            elif line.startswith("DR ") and len(line.strip().split('   ')) > 1:
                # Parse DR (Database Cross-Reference) lines for GO annotations
                items = line.strip().split('   ')[1].split('; ')
                if items[0] == 'GO':
                    # Format: GO; GO:0003674; F:DNA binding; EXP:Inferred from Experiment
                    if len(items) >= 4:
                        evidence_code = items[3].split(':')[0]
                        if evidence_code in EXP_CODES:
                            curr_has_exp_go = True
            elif line.startswith("OX "):
                m = ox_tax_re.search(line)
                if m: curr_taxid = m.group(1)
            elif line.startswith("OS ") and want_names:
                # OS   Homo sapiens (Human).
                name = line[5:].strip()
                # Trim trailing dot and extras split across lines until a line not starting with OS
                name = name.rstrip(".")
                if curr_os is None:
                    curr_os = name
                else:
                    curr_os += " " + name
            elif line.startswith("//"):
                if curr_has_exp_go and curr_taxid:
                    tax_counts[curr_taxid] = tax_counts.get(curr_taxid, 0) + 1
                    if want_names and curr_os:
                        # keep first seen OS as canonical
                        tax_names.setdefault(curr_taxid, curr_os)
                # reset for next entry
                curr_has_exp_go = False
                curr_taxid = None
                curr_os = None
            # else: ignore other lines
    return tax_counts, tax_names

def write_json(path: str, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)

def load_taxid_file(path: str) -> set:
    s = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"): continue
            line = line.split(" ")[0]
            if not line.isdigit():
                print(f"[WARN] Non-numeric taxid ignored: {line}")
                continue
            s.add(line)
    return s

def have_datasets_cli() -> bool:
    try:
        subprocess.run(["datasets", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except Exception:
        return False

def download_ncbi_proteins_for_taxon(taxid: str, out_dir: str, refseq_only: bool, gff: bool = False, is_reference_taxon: bool = False, retries: int = 2) -> Tuple[str, bool, str]:
    """
    Use NCBI Datasets CLI to download genome protein FASTA for a taxon.
    Returns (taxid, success, message)
    """
    tax_dir = os.path.join(out_dir, taxid)
    protein_dir = os.path.join(tax_dir, "protein")
    os.makedirs(protein_dir, exist_ok=True)

    # If already has .faa files, skip
    existing = list(pathlib.Path(protein_dir).rglob("*.faa"))
    if existing:
        # Only skip if it's not exactly one legacy protein.faa that we still need to rename
        if not (len(existing) == 1 and existing[0].name == "protein.faa"):
            return (taxid, True, f"skip (found {len(existing)} .faa files)")

    zip_path = os.path.join(tax_dir, f"{taxid}.zip")
    cmd = ["datasets", "download", "genome", "taxon", taxid, "--include", "protein", "--filename", zip_path, "--no-progressbar"]
    if gff:
        cmd += ["--include", "gff3"]
    if refseq_only:
        cmd += ["--assembly-source", "refseq"]
    # use latest only
    cmd += ["--dehydrated"]  # faster manifest; then rehydrate proteins only
    if is_reference_taxon:
        cmd += ["--reference"]
    try:
        # Step 1: download dehydrated package
        for attempt in range(retries+1):
            rc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if rc.returncode == 0:
                break
            time.sleep(2 * attempt + 1)
        else:
            return (taxid, False, f"datasets download failed: {rc.stderr.decode(errors='ignore')[:200]}")

        # Step 2: rehydrate proteins
        # datasets rehydrate --directory <tax_dir> --include protein
        # Unzip first
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tax_dir)
        # Now rehydrate
        rc2 = subprocess.run(["datasets", "rehydrate", "--directory", tax_dir],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if rc2.returncode != 0:
            return (taxid, False, f"rehydrate failed: {rc2.stderr.decode(errors='ignore')[:200]}")

        # Move all *.faa files into protein_dir, renaming to accession.faa (e.g., GCF_XXXX.faa)
        for p in pathlib.Path(tax_dir).rglob("protein.faa"):  # change download protein.faa to GCF_XXXX.faa
            rel = p.relative_to(tax_dir)
            parts = rel.parts
            accession = None
            if "data" in parts:
                try:
                    accession = parts[parts.index("data") + 1]
                except Exception:
                    accession = None
            if not accession:
                accession = p.parent.name  # fallback
            dst = pathlib.Path(protein_dir) / f"{accession}.faa"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(dst))

        # Move all *.gff files into tax_dir, renaming to accession.gff (e.g., GCF_XXXX.gff)
        for p in pathlib.Path(tax_dir).rglob("*.gff"):
            rel = p.relative_to(tax_dir)
            parts = rel.parts
            accession = None
            if "data" in parts:
                try:
                    accession = parts[parts.index("data") + 1]
                except Exception:
                    accession = None
            if not accession:
                accession = p.parent.name  # fallback
            dst = pathlib.Path(tax_dir) / f"{accession}.gff"
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(dst))

        return (taxid, True, "ok")
    finally:
        # cleanup zip & stray dirs (keep protein_dir)
        if os.path.exists(zip_path):
            try: os.remove(zip_path)
            except: pass
        # remove everything except protein_dir
        for child in pathlib.Path(tax_dir).iterdir():
            if child.name == "protein": continue
            try:
                if child.is_file() and not child.name.endswith(".gff"): child.unlink()
                else: shutil.rmtree(child)
            except Exception:
                pass

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="Output dataset directory")
    ap.add_argument("--max-taxa", type=int, default=None)
    ap.add_argument("--taxid-include", type=str, default=None)
    ap.add_argument("--taxid-exclude", type=str, default=None)
    ap.add_argument("--taxids-file", type=str, default=None, help="File containing taxids to process (one per line). When provided, skips Swiss-Prot processing and only processes these taxa.")
    ap.add_argument("--is-reference-taxon", action="store_true", help="Whether the taxids are reference taxons (e.g., human, mouse, rat, etc.). When provided, uses the --reference flag when downloading NCBI proteomes.")
    ap.add_argument("--gff", action="store_true", help="Whether to download GFF3 files along with the protein FASTA files.")
    ap.add_argument("--refseq-only", action="store_true")
    ap.add_argument("--organism-report", action="store_true")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--skip-ncbi", action="store_true")
    ap.add_argument("--swissprot-url", type=str, default=SWISSPROT_DEFAULT_URL)
    ap.add_argument("--swissprot-strict", action="store_true", help="Force re-download of Swiss-Prot flatfile even if it exists")
    args = ap.parse_args()

    out_dir = os.path.abspath(args.out)
    data_dir = os.path.join(out_dir, "data")
    meta_dir = os.path.join(out_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    today = datetime.date.today().isoformat()
    
    if args.taxids_file:
        # Skip Swiss-Prot processing and use taxids from file
        print(f"[TaxIDs] Loading taxids from file: {args.taxids_file}")
        taxa = load_taxid_file(args.taxids_file)
        print(f"[TaxIDs] Loaded {len(taxa)} taxids from file")
        
        # Create dummy counts for compatibility (all set to 1)
        counts = {taxid: 1 for taxid in taxa}
        names = {}  # No organism names available when using taxids file
        
        # Sort taxa alphabetically since we don't have protein counts
        taxa_sorted = sorted(taxa)
        if args.max_taxa is not None:
            taxa_sorted = taxa_sorted[:args.max_taxa]

    else:
        # Original Swiss-Prot processing
        # 1) Download Swiss-Prot flatfile (skip if exists unless --strict)
        sp_dir = os.path.join(data_dir, "uniprot")
        os.makedirs(sp_dir, exist_ok=True)
        sp_path = os.path.join(sp_dir, "uniprot_sprot.dat.gz")
        
        if os.path.exists(sp_path) and not args.swissprot_strict:
            print(f"[Swiss-Prot] File already exists: {sp_path}")
            print("[Swiss-Prot] Skipping download (use --swissprot-strict to force re-download)")
        else:
            print(f"[Swiss-Prot] Downloading: {args.swissprot_url}")
            download(args.swissprot_url, sp_path)
        
            sp_checksum = sha256_of_file(sp_path)
            sp_meta = {
                "source": args.swissprot_url,
                "download_date": today,
                "filename": os.path.basename(sp_path),
                "checksum_sha256": sp_checksum
            }
            write_json(os.path.join(meta_dir, f"uniprot_sprot_{today}.json"), sp_meta)

        # 2) Parse for experimental GO annotations & taxa
        print("[Swiss-Prot] Parsing experimental GO annotations to collect NCBI TaxIDs...")
        counts, names = parse_uniprot_dat_for_taxa(sp_path, want_names=args.organism_report)
        print(f"[Swiss-Prot] Found {sum(counts.values())} proteins with experimental GO annotations across {len(counts)} taxa.")

        # 3) Apply include/exclude filters
        taxa = set(counts.keys())
        if args.taxid_include:
            include = load_taxid_file(args.taxid_include)
            taxa = taxa.intersection(include)
        if args.taxid_exclude:
            exclude = load_taxid_file(args.taxid_exclude)
            taxa = taxa.difference(exclude)

        # Sort taxa by number of proteins desc
        taxa_sorted = sorted(taxa, key=lambda t: -counts[t])
        if args.max_taxa is not None:
            taxa_sorted = taxa_sorted[:args.max_taxa]

        # Save taxon list
        tax_list_path = os.path.join(out_dir, f"taxa_from_swissprot_exp_go_{today}.tsv")
        with open(tax_list_path, "w", encoding="utf-8") as f:
            header = ["taxid", "n_exp_go_proteins"]
            if args.organism_report: header.append("organism")
            f.write("\t".join(header) + "\n")
            for t in taxa_sorted:
                row = [t, str(counts[t])]
                if args.organism_report:
                    row.append(names.get(t, ""))
                f.write("\t".join(row) + "\n")
        print(f"[Swiss-Prot] Wrote taxon list: {tax_list_path}")

    # 4) Download NCBI proteomes for these taxa
    if args.skip_ncbi:
        print("[NCBI] Skipping NCBI downloads as requested (--skip-ncbi).")
        return

    if not have_datasets_cli():
        print("[ERROR] NCBI Datasets CLI (`datasets`) not found in PATH.", file=sys.stderr)
        print("Install instructions: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/", file=sys.stderr)
        sys.exit(2)

    ncbi_dir = os.path.join(data_dir, "ncbi_proteins")
    os.makedirs(ncbi_dir, exist_ok=True)

    print(f"[NCBI] Downloading proteomes for {len(taxa_sorted)} taxa using NCBI Datasets (threads={args.threads})...")
    results = []
    with cf.ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futs = [ex.submit(download_ncbi_proteins_for_taxon, t, ncbi_dir, args.refseq_only, args.gff, is_reference_taxon=args.is_reference_taxon) for t in taxa_sorted]
        for fut in cf.as_completed(futs):
            taxid, ok, msg = fut.result()
            results.append((taxid, ok, msg))
            status = "OK" if ok else "FAIL"
            print(f"[NCBI] {taxid}: {status} - {msg}")

    # 5) Write a run report & checksums
    report = {
        "run_date": today,
        "n_taxa_requested": len(taxa_sorted),
        "n_taxa_ok": sum(1 for _,ok,_ in results if ok),
        "n_taxa_fail": sum(1 for _,ok,_ in results if not ok),
        "refseq_only": args.refseq_only,
        "swissprot_source": args.swissprot_url
    }
    write_json(os.path.join(out_dir, f"run_report_{today}.json"), report)

    # Checksums for NCBI .faa files
    checksum_path = os.path.join(out_dir, f"ncbi_proteins_sha256_{today}.tsv")
    with open(checksum_path, "w", encoding="utf-8") as f:
        f.write("sha256\tpath\n")
        for p in pathlib.Path(ncbi_dir).rglob("*.faa"):
            s = sha256_of_file(str(p))
            rel = os.path.relpath(str(p), out_dir)
            f.write(f"{s}\t{rel}\n")
    print(f"[DONE] Report: {os.path.join(out_dir, f'run_report_{today}.json')}")
    print(f"[DONE] Checksums: {checksum_path}")
    print("[DONE] Dataset is ready.")
    
if __name__ == "__main__":
    main()
