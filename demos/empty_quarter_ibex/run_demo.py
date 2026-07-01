#!/usr/bin/env python3
"""Run GenoAdjust on a real Empty Quarter genome annotation.

The demo uses PGAP's GO-bearing CDS annotations for rh04 / 63_rh04 and converts
them to the repository's tabular prediction format. PGAP does not emit
probabilities in this GFF, so each observed GO annotation is encoded with the
same high score and thresholded in the usual way by the optimizer. By default,
inputs are fetched from the public bio2vec.net reproducibility bundle; use
--source ibex to refresh from the original IBEX paths.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote
from urllib.request import urlopen


REPO_ROOT = Path(__file__).resolve().parents[2]
DEMO_ROOT = Path(__file__).resolve().parent

SAMPLE_ID = "rh04"
GENOME_ID = "63_rh04"
TAXON_ID = "96241"
TAXON_NAME = "Bacillus spizizenii"
THRESHOLD = 0.5
DEFAULT_SCORE = 0.9

REMOTE_ROOT = (
    "/ibex/scratch/projects/c2014/EmptyQuarter_Data/"
    "cultures/assemblies/assemblies_tiannyu_2025/pgap_annotations/63_rh04"
)
REMOTE_GFF = f"{REMOTE_ROOT}/annot.gff"
REMOTE_ANI_REPORT = f"{REMOTE_ROOT}/ani-tax-report.txt"

PUBLIC_BASE_URL = (
    "https://bio2vec.net/data/genoadjust/"
    "rh04_63_rh04_bacillus_spizizenii_genoadjust_demo"
)
PUBLIC_GFF_URL = f"{PUBLIC_BASE_URL}/pgap_annot.gff.gz"
PUBLIC_ANI_REPORT_URL = f"{PUBLIC_BASE_URL}/ani-tax-report.txt"

BACILLUS_LINEAGE = [
    ("NCBITaxon_96241", "NCBITaxon_653685"),
    ("NCBITaxon_653685", "NCBITaxon_1386"),
    ("NCBITaxon_1386", "NCBITaxon_186817"),
    ("NCBITaxon_186817", "NCBITaxon_1385"),
    ("NCBITaxon_1385", "NCBITaxon_91061"),
    ("NCBITaxon_91061", "NCBITaxon_1239"),
    ("NCBITaxon_1239", "NCBITaxon_1783272"),
    ("NCBITaxon_1783272", "NCBITaxon_2"),
    ("NCBITaxon_2", "NCBITaxon_131567"),
    ("NCBITaxon_131567", "NCBITaxon_1"),
]

GO_RE = re.compile(r"GO:\d{7}")
GO_LABEL_RE = re.compile(r"([^|]+)\|(\d{7})\|\|[A-Z]+")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch rh04 PGAP annotations from IBEX and run GenoAdjust."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEMO_ROOT / "output",
        help="Directory for fetched inputs, converted predictions, logs, and adjusted outputs.",
    )
    parser.add_argument(
        "--source",
        choices=("public", "ibex"),
        default="public",
        help="Fetch demo inputs from bio2vec.net/data by default, or refresh them from IBEX.",
    )
    parser.add_argument(
        "--base-url",
        default=PUBLIC_BASE_URL,
        help="Public base URL for the reproducibility bundle when --source=public.",
    )
    parser.add_argument(
        "--remote-host",
        default="ibex",
        help="SSH host alias used by scp to fetch the IBEX annotation files.",
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Reuse output/raw/annot.gff and output/raw/ani-tax-report.txt instead of fetching them.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=2,
        help="Candidate proteins per singleton complex for the optimized Stage 2 solver.",
    )
    return parser.parse_args()


def parse_gff_attributes(field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for part in field.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def download_url(url: str, path: Path) -> None:
    with urlopen(url, timeout=60) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def fetch_inputs(
    output_dir: Path,
    source: str,
    base_url: str,
    remote_host: str,
    skip_fetch: bool,
) -> tuple[Path, Path, str]:
    raw_dir = output_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    gff_path = raw_dir / "annot.gff"
    ani_path = raw_dir / "ani-tax-report.txt"
    if skip_fetch:
        missing = [str(p) for p in (gff_path, ani_path) if not p.exists()]
        if missing:
            raise SystemExit(f"--skip-fetch requested but files are missing: {', '.join(missing)}")
        return gff_path, ani_path, "local cached files"

    if source == "public":
        base_url = base_url.rstrip("/")
        gff_gz_path = raw_dir / "pgap_annot.gff.gz"
        download_url(f"{base_url}/pgap_annot.gff.gz", gff_gz_path)
        download_url(f"{base_url}/ani-tax-report.txt", ani_path)
        with gzip.open(gff_gz_path, "rb") as src, gff_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        return gff_path, ani_path, f"{base_url}/pgap_annot.gff.gz"

    for remote_path, local_path in (
        (REMOTE_GFF, gff_path),
        (REMOTE_ANI_REPORT, ani_path),
    ):
        subprocess.run(
            ["scp", "-q", f"{remote_host}:{remote_path}", str(local_path)],
            check=True,
        )
    return gff_path, ani_path, REMOTE_GFF


def convert_gff_to_predictions(
    gff_path: Path,
    predictions_path: Path,
    metadata_path: Path,
    score: float = DEFAULT_SCORE,
) -> dict[str, object]:
    protein_terms: dict[str, set[str]] = defaultdict(set)
    term_labels: dict[str, str] = {}
    cds_records = 0
    cds_with_go = 0

    with gff_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            cds_records += 1
            attrs = parse_gff_attributes(fields[8])
            terms = set(GO_RE.findall(attrs.get("Ontology_term", "")))
            if not terms:
                continue
            cds_with_go += 1
            protein_id = attrs.get("locus_tag") or attrs.get("protein_id") or attrs.get("ID")
            if not protein_id:
                continue
            protein_terms[protein_id].update(terms)
            for key in ("go_function", "go_process", "go_component"):
                value = attrs.get(key, "")
                for match in GO_LABEL_RE.finditer(value):
                    label = match.group(1).lstrip(", ").strip()
                    term = f"GO:{match.group(2)}"
                    if label:
                        term_labels.setdefault(term, label)

    predictions_path.parent.mkdir(parents=True, exist_ok=True)
    with predictions_path.open("w") as handle:
        for protein_id in sorted(protein_terms):
            encoded_terms = [
                f"{term}|{score:.6f}" for term in sorted(protein_terms[protein_id])
            ]
            handle.write("\t".join([protein_id, *encoded_terms]) + "\n")

    annotation_count = sum(len(terms) for terms in protein_terms.values())
    metadata = {
        "sample_id": SAMPLE_ID,
        "genome_id": GENOME_ID,
        "taxon_id": TAXON_ID,
        "taxon_name": TAXON_NAME,
        "remote_gff": REMOTE_GFF,
        "cds_records": cds_records,
        "cds_with_go": cds_with_go,
        "proteins_with_go": len(protein_terms),
        "distinct_go_terms": len({term for terms in protein_terms.values() for term in terms}),
        "prediction_annotations": annotation_count,
        "score": score,
        "threshold": THRESHOLD,
        "term_labels": term_labels,
    }
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def write_demo_lineage(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for child, parent in BACILLUS_LINEAGE:
            handle.write(f"{child}\t{parent}\n")


def run_logged(command: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )
    log_path.write_text(
        "$ " + shlex.join(command) + "\n\n" + proc.stdout + proc.stderr
    )
    if proc.returncode != 0:
        sys.stderr.write(log_path.read_text())
        raise SystemExit(proc.returncode)


def load_predictions(path: Path) -> dict[str, dict[str, float]]:
    predictions: dict[str, dict[str, float]] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            protein_id = fields[0]
            predictions[protein_id] = {}
            for term_score in fields[1:]:
                term, score = term_score.split("|", 1)
                predictions[protein_id][term] = float(score)
    return predictions


def above_threshold(predictions: dict[str, dict[str, float]]) -> set[tuple[str, str]]:
    return {
        (protein_id, term)
        for protein_id, terms in predictions.items()
        for term, score in terms.items()
        if score > THRESHOLD
    }


def parse_ani_report(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.startswith("Submitted organism:"):
            result["submitted"] = line.split(":", 1)[1].strip()
        elif line.startswith("Best match:"):
            result["best_match"] = line.split(":", 1)[1].strip()
        elif line.startswith("Status:"):
            result["status"] = line.split(":", 1)[1].strip()
        elif line.startswith("Confidence:"):
            result["confidence"] = line.split(":", 1)[1].strip()
    return result


def format_term_counts(
    pairs: set[tuple[str, str]],
    labels: dict[str, str],
    limit: int,
) -> list[str]:
    counts = Counter(term for _, term in pairs)
    lines = []
    for term, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:limit]:
        label = labels.get(term, "unlabeled in PGAP fields")
        lines.append(f"  {term} ({label}): {count}")
    return lines


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    raw_gff, ani_report, source_gff = fetch_inputs(
        output_dir,
        source=args.source,
        base_url=args.base_url,
        remote_host=args.remote_host,
        skip_fetch=args.skip_fetch,
    )

    predictions_dir = output_dir / "predictions"
    optimized_dir = output_dir / "optimized"
    constraints_dir = output_dir / "constraints"
    logs_dir = output_dir / "logs"
    prediction_file = predictions_dir / f"predictions_fold_00_taxon_{TAXON_ID}.tsv"
    taxon_output = optimized_dir / f"taxon_fold_00_taxon_{TAXON_ID}.tsv"
    final_output = optimized_dir / f"optimized_fold_00_taxon_{TAXON_ID}.tsv"
    metadata_path = output_dir / "metadata.json"
    lineage_path = constraints_dir / "ncbitaxon_hierarchy.tsv"

    optimized_dir.mkdir(parents=True, exist_ok=True)
    metadata = convert_gff_to_predictions(raw_gff, prediction_file, metadata_path)
    write_demo_lineage(lineage_path)

    taxon_command = [
        sys.executable,
        "-u",
        "taxon_consistency/adjust_ortools.py",
        "--predictions",
        str(prediction_file),
        "--constraints",
        str(REPO_ROOT / "data/go_taxon_constraints_extracted_obo.tsv"),
        "--go-hierarchy",
        str(REPO_ROOT / "data/go_hierarchy.tsv"),
        "--taxon-hierarchy",
        str(REPO_ROOT / "data/taxon_hierarchy.tsv"),
        "--ncbitaxon-hierarchy",
        str(lineage_path),
        "--output",
        str(taxon_output),
        "--threshold",
        str(THRESHOLD),
        "--taxon-id",
        TAXON_ID,
    ]
    run_logged(taxon_command, logs_dir / "stage1_taxon_consistency.log")

    complex_command = [
        sys.executable,
        "-u",
        "complex_coherence/adjust_ortools.py",
        "--predictions",
        str(taxon_output),
        "--complexes",
        str(REPO_ROOT / "complex_coherence/protein_complexes.tsv"),
        "--go_hierarchy",
        str(REPO_ROOT / "data/go_hierarchy.tsv"),
        "--output",
        str(final_output),
        "--threshold",
        str(THRESHOLD),
        "--optimized",
        "--top_k",
        str(args.top_k),
    ]
    run_logged(complex_command, logs_dir / "stage2_complex_coherence.log")

    before = above_threshold(load_predictions(prediction_file))
    after_taxon = above_threshold(load_predictions(taxon_output))
    after_complex = above_threshold(load_predictions(final_output))
    taxon_demotions = before - after_taxon
    taxon_promotions = after_taxon - before
    complex_demotions = after_taxon - after_complex
    complex_promotions = after_complex - after_taxon
    labels = metadata["term_labels"]
    assert isinstance(labels, dict)
    ani = parse_ani_report(ani_report)

    print("Source genome")
    print(f"  sample: {SAMPLE_ID} / {GENOME_ID}")
    print(f"  GFF source: {source_gff}")
    print(f"  ANI best match: {ani.get('best_match', 'not reported')}")
    print(f"  ANI status: {ani.get('status', 'not reported')} ({ani.get('confidence', 'not reported')})")
    print()
    print("Converted PGAP GO predictions")
    print(f"  CDS records: {metadata['cds_records']}")
    print(f"  CDS records with GO terms: {metadata['cds_with_go']}")
    print(f"  proteins with GO predictions: {metadata['proteins_with_go']}")
    print(f"  distinct GO terms: {metadata['distinct_go_terms']}")
    print(f"  above-threshold annotations before adjustment: {len(before)}")
    print()
    print("Stage 1: taxon consistency")
    print(f"  above-threshold annotations after Stage 1: {len(after_taxon)}")
    print(f"  demotions: {len(taxon_demotions)}")
    print(f"  promotions: {len(taxon_promotions)}")
    for line in format_term_counts(taxon_demotions, labels, limit=8):
        print(line)
    print()
    print("Stage 2: complex coherence")
    print(f"  top_k: {args.top_k}")
    print(f"  above-threshold annotations after Stage 2: {len(after_complex)}")
    print(f"  demotions: {len(complex_demotions)}")
    print(f"  promotions: {len(complex_promotions)}")
    for line in format_term_counts(complex_demotions, labels, limit=8):
        print(line)
    print()
    print("Generated files")
    print(f"  predictions: {prediction_file.relative_to(REPO_ROOT)}")
    print(f"  taxon-adjusted: {taxon_output.relative_to(REPO_ROOT)}")
    print(f"  final adjusted: {final_output.relative_to(REPO_ROOT)}")
    print(f"  raw solver logs: {logs_dir.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
