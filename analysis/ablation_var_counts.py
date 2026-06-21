#!/usr/bin/env python3
"""
Materialize the Stage-2 scalability ablation (paper Section 3.6) to a CSV.

For each of the three benchmark proteomes (P. abyssi 272844 small, 759272
medium, R. norvegicus 10116 large) this records, per heuristic configuration,
the number of Boolean variables and constraints the model would contain
(dry-run counts, so the full Naive configuration does not OOM) and the wall-clock
solve time of the finishing +H1+H2+H3 configuration. This is the committed
artefact behind the headline reduction (164M / 216M variables/constraints down
to ~30,000 each, a 5,380x reduction, solved in ~13 s).

Usage:
    python analysis/ablation_var_counts.py \
        --predictions_dir ${DATA_DIR}/swissprot_proteomes_folds/seq_sim_timeset_results/cc/predictions \
        --output analysis/revision_results/ablation_var_counts.csv
"""

import argparse
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from complex_coherence.heuristic_ablation import run_ablation  # noqa: E402

BENCHMARKS = [
    ("272844", "Pyrococcus abyssi (small)"),
    ("759272", "medium"),
    ("10116", "Rattus norvegicus (large)"),
]


def find_predictions(predictions_dir, taxon_id):
    hits = sorted(glob.glob(os.path.join(
        predictions_dir, f"predictions_fold_*_taxon_{taxon_id}.tsv")))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions_dir", required=True)
    ap.add_argument("--complexes", default="complex_coherence/protein_complexes.tsv")
    ap.add_argument("--go_hierarchy", default="data/go_hierarchy_cc.tsv")
    ap.add_argument("--threshold", type=float, default=0.3)
    ap.add_argument("--top_k", type=int, default=50)
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--output", default="analysis/revision_results/ablation_var_counts.csv")
    ap.add_argument("--solve", action="store_true",
                    help="Also time the +H1+H2+H3 solve for every proteome (slow).")
    ap.add_argument("--solve_taxon", default="10116",
                    help="Always time the +H1+H2+H3 solve for this taxon (the large proteome).")
    args = ap.parse_args()

    rows = []
    for taxon_id, organism in BENCHMARKS:
        pred = find_predictions(args.predictions_dir, taxon_id)
        if not pred:
            print(f"[WARN] no predictions for taxon {taxon_id}, skipping")
            continue

        # Dry-run: variable/constraint counts for every configuration (no model
        # is built, so the Naive full-variable config cannot OOM).
        dry = run_ablation(pred, args.complexes, args.go_hierarchy, taxon_id,
                           args.threshold, args.top_k, args.timeout, dry_run=True)
        # Real solve of the finishing configuration to record solve_s. The
        # preprocessing is expensive on the large proteome, so only solve where
        # requested (default: the large proteome, which carries the headline).
        solve_s, status = {}, {}
        if args.solve or taxon_id == args.solve_taxon:
            solved = run_ablation(pred, args.complexes, args.go_hierarchy, taxon_id,
                                  args.threshold, args.top_k, args.timeout,
                                  dry_run=False, configs={"+H1+H2+H3"})
            solve_s = {r.config: r.solve_s for r in solved}
            status = {r.config: r.status for r in solved}

        for r in dry:
            rows.append({
                "taxon_id": taxon_id,
                "organism": organism,
                "config": r.config,
                "n_vars": r.n_vars,
                "n_constraints": r.n_constraints,
                "setup_s": f"{r.setup_s:.3f}",
                "solve_s": f"{solve_s.get(r.config, ''):.3f}" if r.config in solve_s else "",
                "solved_status": status.get(r.config, ""),
            })

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "taxon_id", "organism", "config", "n_vars", "n_constraints",
            "setup_s", "solve_s", "solved_status"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows to {args.output}")
    # Echo the headline reduction for the large proteome.
    for r in rows:
        if r["taxon_id"] == "10116" and r["config"] in ("Naive", "+H1+H2+H3"):
            print(f"  {r['organism']} {r['config']}: "
                  f"vars={r['n_vars']} constraints={r['n_constraints']} solve_s={r['solve_s']}")


if __name__ == "__main__":
    main()
