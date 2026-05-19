"""
Match up a proteome's Uniprot protein identifiers to EMBL ids and their locations in the genome.

Pipeline:
  1. Parse .dat file  -> protein_id -> (entry_name, is_swissprot)
  2. Parse genome file -> protein_id -> (component, start, end, is_complement)
     Coordinates are made flat by concatenating all components in order of appearance.
  3. Join on protein_id and write a TSV compatible with metagenomics_plot.py.

Output TSV columns: start, end, complement, entry_name, is_swissprot, component
  - 'start' and 'end' are 1-based, inclusive, in concatenated genome space (offset by the component's size).
  - Columns 0-3 match the format expected by metagenomics_plot.py's load_cds_positions().
"""

import re
import os
import argparse
import sys
from collections import OrderedDict
from typing import Dict, Tuple, List
import gzip


# ---------------------------------------------------------------------------
# Step 1: Parse the UniProt .dat file
# ---------------------------------------------------------------------------

def parse_dat_file(dat_path: str) -> Dict[str, Tuple[str, bool]]:
    """
    Parse a UniProt flat-file (.dat or .dat.gz) and return a mapping:
        protein_id -> (entry_name, is_swissprot)

    - entry_name is the first token of the ID line (e.g. 'ATG9_DROME').
    - protein_id is the EMBL protein accession (e.g. 'AGB93555.1') from
      'DR   EMBL; ...' lines where the data-class is 'Genomic_DNA'.
    - is_swissprot is True when the ID line says 'Reviewed'.
    - When multiple UniProt entries claim the same protein_id the first
      encountered wins (a warning is printed for the collision).
    """
    protein_id_to_entry: Dict[str, Tuple[str, bool]] = {}

    embl_dr_re = re.compile(
        r'^DR\s+EMBL;\s+\S+;\s+(\S+);\s+[^;]+;\s+(\S+)\.'
    )

    current_entry_name: str = ""
    current_is_swissprot: bool = False
    current_protein_ids: List[str] = []

    def flush():
        for pid in current_protein_ids:
            if pid == "-":
                continue
            if pid in protein_id_to_entry:
                existing = protein_id_to_entry[pid][0]
                if existing != current_entry_name:
                    print(
                        f"[WARN] protein_id {pid} already mapped to "
                        f"{existing}; ignoring {current_entry_name}"
                    )
            else:
                protein_id_to_entry[pid] = (current_entry_name, current_is_swissprot)

    opener = gzip.open if dat_path.endswith(".gz") else open
    with opener(dat_path, "rt") as fh:
        for line in fh:
            if line.startswith("ID   "):
                # 'ID   ATG9_DROME              Reviewed;         852 AA.'
                # 'ID   A0A0B4KF86_DROME        Unreviewed;       852 AA.'
                tokens = line.split()
                current_entry_name = tokens[1] if len(tokens) >= 2 else ""
                current_is_swissprot = "Reviewed;" in line and "Unreviewed;" not in line
                current_protein_ids = []

            elif line.startswith("DR   EMBL;"):
                m = embl_dr_re.match(line)
                if m:
                    pid = m.group(1).strip()
                    data_class = m.group(2).strip()
                    if data_class == "Genomic_DNA":
                        current_protein_ids.append(pid)

            elif line.startswith("//"):
                if current_entry_name:
                    flush()
                current_entry_name = ""
                current_is_swissprot = False
                current_protein_ids = []

    print(
        f"[INFO] .dat: mapped {len(protein_id_to_entry)} EMBL protein_ids "
        f"to UniProt entry names"
    )
    return protein_id_to_entry


# ---------------------------------------------------------------------------
# Step 2: Parse the genome EMBL flat file
# ---------------------------------------------------------------------------

_COORD_RE = re.compile(r"\d+")


def _parse_location(loc_str: str) -> Tuple[int, int, bool]:
    """
    Given a (possibly multi-exon) EMBL location string return
    (min_start, max_end, is_complement) where coordinates are 1-based.

    Handles:
      75280..76098
      complement(115713..116117)
      join(67625..67762,67892..68023,68085..70549,70607..70895)
      complement(join(11215..11344,11410..11518,...))
    """
    is_complement = loc_str.startswith("complement(")
    nums = [int(n) for n in _COORD_RE.findall(loc_str)]
    if not nums:
        raise ValueError(f"No coordinates found in location: {loc_str!r}")
    return min(nums), max(nums), is_complement


def parse_genome_file(
    genome_path: str,
) -> Tuple[Dict[str, Tuple[str, int, int, bool]], OrderedDict]:
    """
    Parse an EMBL genome flat file (potentially containing multiple components).

    Returns:
        protein_id_to_cds: protein_id -> (component_accession, start, end, is_complement)
            Coordinates are LOCAL to the component (1-based, inclusive).
        component_sizes: OrderedDict[component_accession -> size_bp]
            Ordered by first appearance in the file.
    """
    protein_id_re = re.compile(r'/protein_id="([^"]+)"')
    id_line_re = re.compile(
        r'^ID\s+(\S+?)(?:;|\s).*?;\s*(\d+)\s+BP\.'
    )

    protein_id_to_cds: Dict[str, Tuple[str, int, int, bool]] = {}
    component_sizes: OrderedDict = OrderedDict()

    current_component: str = ""
    current_size: int = 0

    # CDS block state
    in_cds: bool = False
    loc_lines: List[str] = ""   # will be a list
    loc_lines = []
    loc_done: bool = False      # True once we've moved past the location lines
    current_loc: str = ""
    current_pid: str = ""

    def flush_cds():
        nonlocal current_loc, current_pid
        if not current_loc:
            return
        try:
            start, end, is_comp = _parse_location(current_loc)
        except ValueError as e:
            print(f"[WARN] Could not parse location {current_loc!r}: {e}")
            current_loc = ""
            current_pid = ""
            return
        if current_pid:
            if current_pid in protein_id_to_cds:
                pass  # keep first occurrence
            else:
                protein_id_to_cds[current_pid] = (
                    current_component, start, end, is_comp
                )
        current_loc = ""
        current_pid = ""

    new_feature_re = re.compile(r"^FT   \S")
    cds_start_re = re.compile(r"^FT   CDS\s+(.+)$")
    loc_cont_re = re.compile(r"^FT\s{19,}([^/].*)$")   # continuation, not a qualifier

    with open(genome_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # --- Component header ---
            if line.startswith("ID   "):
                # Flush any pending CDS from previous component
                if in_cds:
                    flush_cds()
                    in_cds = False

                m = id_line_re.match(line)
                if m:
                    current_component = m.group(1)
                    current_size = int(m.group(2))
                    if current_component not in component_sizes:
                        component_sizes[current_component] = current_size
                else:
                    print(f"[WARN] Could not parse ID line: {line!r}")
                    current_component = ""
                    current_size = 0
                continue

            # --- Feature table ---
            if line.startswith("FT   CDS"):
                # Start of a new CDS feature
                if in_cds:
                    flush_cds()
                m = cds_start_re.match(line)
                loc_raw = m.group(1).strip() if m else ""
                in_cds = True
                loc_done = False
                current_loc = loc_raw
                current_pid = ""
                continue

            if in_cds:
                # Check for a new (non-continuation) feature line
                if new_feature_re.match(line) and not line.startswith("FT   CDS"):
                    flush_cds()
                    in_cds = False
                    continue

                # Location continuation (before any qualifier)
                if not loc_done:
                    cont = loc_cont_re.match(line)
                    if cont:
                        current_loc += cont.group(1).strip()
                        continue
                    else:
                        # We've hit a qualifier line (/...) or end of FT block
                        loc_done = True

                # Qualifier lines
                pid_m = protein_id_re.search(line)
                if pid_m:
                    current_pid = pid_m.group(1)

        # End of file
        if in_cds:
            flush_cds()

    print(
        f"[INFO] genome: found {len(component_sizes)} components, "
        f"{len(protein_id_to_cds)} CDS entries with protein_id"
    )
    return protein_id_to_cds, component_sizes


# ---------------------------------------------------------------------------
# Step 3: Match and write output
# ---------------------------------------------------------------------------

def build_output(
    protein_id_to_entry: Dict[str, Tuple[str, bool]],
    protein_id_to_cds: Dict[str, Tuple[str, int, int, bool]],
    component_sizes: OrderedDict,
    output_path: str,
) -> int:
    """
    Join the two mappings on protein_id, apply per-component coordinate offsets,
    and write the output TSV.

    Returns the number of rows written.
    """
    # Compute cumulative offsets for each component (order = appearance in genome file)
    offsets: Dict[str, int] = {}
    cumulative = 0
    for comp, size in component_sizes.items():
        offsets[comp] = cumulative
        cumulative += size

    rows_written = 0
    duplicates_skipped = 0
    missing_in_genome = 0
    missing_in_dat = 0

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    seen: set = set()

    with open(output_path, "w") as out:
        out.write("start\tend\tcomplement\tentry_name\tis_swissprot\tcomponent\n")

        for pid, (entry_name, is_sp) in sorted(protein_id_to_entry.items()):
            cds = protein_id_to_cds.get(pid)
            if cds is None:
                missing_in_genome += 1
                continue
            comp, start, end, is_comp = cds
            offset = offsets.get(comp, 0)
            row_key = (start + offset, end + offset, is_comp, entry_name, comp)
            if row_key in seen:
                duplicates_skipped += 1
                continue
            seen.add(row_key)
            out.write(
                f"{start + offset}\t{end + offset}\t{is_comp}\t"
                f"{entry_name}\t{is_sp}\t{comp}\n"
            )
            rows_written += 1

        for pid in protein_id_to_cds:
            if pid not in protein_id_to_entry:
                missing_in_dat += 1

    print(f"[INFO] Rows written:              {rows_written}")
    print(f"[INFO] Duplicate rows skipped:    {duplicates_skipped}")
    print(f"[INFO] protein_ids not in genome: {missing_in_genome}")
    print(f"[INFO] protein_ids not in .dat:   {missing_in_dat}")
    return rows_written


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Match UniProt entry names to CDS genomic locations via EMBL protein_ids. "
            "Handles multi-chromosome / multi-component genome files."
        )
    )
    parser.add_argument(
        "--dat-file", required=True,
        help="UniProt flat-file (.dat) for the proteome, e.g. UP000000803_7227.dat"
    )
    parser.add_argument(
        "--genome-file", required=True,
        help="EMBL genome flat-file, e.g. GCA_000001215.4.txt"
    )
    parser.add_argument(
        "--output-file", default="CDS_swissprot.tsv",
        help="Output TSV path (default: CDS_swissprot.tsv)"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Step 1: Parsing .dat file ...")
    print("=" * 60)
    protein_id_to_entry = parse_dat_file(args.dat_file)

    print()
    print("=" * 60)
    print("Step 2: Parsing genome file ...")
    print("=" * 60)
    protein_id_to_cds, component_sizes = parse_genome_file(args.genome_file)

    total_genome_bp = sum(component_sizes.values())
    print(
        f"[INFO] Total concatenated genome size: {total_genome_bp:,} bp "
        f"across {len(component_sizes)} components"
    )

    print()
    print("=" * 60)
    print(f"Step 3: Matching and writing to {args.output_file} ...")
    print("=" * 60)
    n = build_output(
        protein_id_to_entry, protein_id_to_cds, component_sizes, args.output_file
    )

    if n == 0:
        print(
            "[ERROR] No rows written. Check that the .dat and genome files "
            "correspond to the same proteome."
        )
        sys.exit(1)

    print()
    print(f"Done. {n} CDS entries written to {args.output_file}")


if __name__ == "__main__":
    main()
