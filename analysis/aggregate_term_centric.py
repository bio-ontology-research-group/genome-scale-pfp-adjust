#!/usr/bin/env python3
"""
aggregate_term_centric.py

Reads run_log.txt produced by analysis/run_term_centric_timeset.sh and all per-term
detail TSVs in term-centric-analysis-v2/ to produce two thesis-ready outputs:

1. flip_precision_table.tsv
   One row per method × ontology × condition.
   Columns: method, ont, condition, n_terms_changed, fps_removed, tps_lost,
            flip_precision, tmax

2. binned_summary_combined.tsv
   All binned summaries stacked with method/ont/cond prefix columns.

Run from: repository root
  python3 analysis/aggregate_term_centric.py
"""

import csv, os, re, sys

OUT_BASE = "term-centric-analysis-v2"
LOG_FILE = os.path.join(OUT_BASE, "run_log.txt")

METHODS  = ["seq_sim", "deepgo", "sprof"]
ONTS     = ["cc", "mf", "bp"]
CONDS    = ["taxons", "without_taxons"]

METHOD_LABEL = {"seq_sim": "Seq-Sim", "deepgo": "DeepGO-SE", "sprof": "SPROF-GO"}
ONT_LABEL    = {"cc": "CC", "mf": "MF", "bp": "BP"}
COND_LABEL   = {"taxons": "taxon-guided", "without_taxons": "taxon-unknown"}


# ── Parse flip-precision from run_log.txt ─────────────────────────────────────

def parse_log(log_path):
    """Return dict (method, ont, cond) -> {fps_removed, tps_lost, flip_precision, tmax}."""
    results = {}
    tag = None
    with open(log_path) as f:
        for line in f:
            line = line.rstrip()
            m = re.match(r'=== (\w+)_(\w+)_(\w+) ===', line)
            if m:
                tag = (m.group(1), m.group(2), m.group(3))
                continue
            if tag and line.startswith("Terms changed"):
                # "Terms changed by solver: N  (FPs removed: X, TPs lost: Y)  Flip precision: Z  ..."
                n_m = re.search(r'Terms changed by solver: (\d+)', line)
                fp_m = re.search(r'FPs removed: (\d+)', line)
                tp_m = re.search(r'TPs lost: (\d+)', line)
                pr_m = re.search(r'Flip precision: ([0-9.]+|nan)', line)
                results[tag] = {
                    'n_changed': int(n_m.group(1)) if n_m else None,
                    'fps_removed': int(fp_m.group(1)) if fp_m else None,
                    'tps_lost': int(tp_m.group(1)) if tp_m else None,
                    'flip_precision': pr_m.group(1) if pr_m else None,
                }
            if tag and line.startswith("F-max threshold"):
                tmax_m = re.search(r'tmax\): ([0-9.]+)', line)
                if tmax_m and tag in results:
                    results[tag]['tmax'] = tmax_m.group(1)
    return results


# ── Write flip-precision table ────────────────────────────────────────────────

def write_flip_precision_table(results, out_path):
    cols = ['method', 'ont', 'condition', 'tmax',
            'n_terms_changed', 'fps_removed', 'tps_lost', 'flip_precision']
    rows = []
    for method in METHODS:
        for ont in ONTS:
            for cond in CONDS:
                tag = (method, ont, cond)
                r = results.get(tag, {})
                rows.append({
                    'method': METHOD_LABEL.get(method, method),
                    'ont': ONT_LABEL.get(ont, ont),
                    'condition': COND_LABEL.get(cond, cond),
                    'tmax': r.get('tmax', ''),
                    'n_terms_changed': r.get('n_changed', ''),
                    'fps_removed': r.get('fps_removed', ''),
                    'tps_lost': r.get('tps_lost', ''),
                    'flip_precision': r.get('flip_precision', ''),
                })
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_path}")


# ── Stack binned summaries ────────────────────────────────────────────────────

def write_binned_combined(out_path):
    all_rows = []
    for method in METHODS:
        for ont in ONTS:
            for cond in CONDS:
                tag = f"{method}_{ont}_{cond}"
                detail_dir = os.path.join(OUT_BASE, f"{tag}_detail")
                binned_file = os.path.join(detail_dir, f"binned_summary_{ont}.tsv")
                if not os.path.exists(binned_file):
                    continue
                with open(binned_file) as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        if not row.get('n_terms'):
                            continue
                        all_rows.append({
                            'method': METHOD_LABEL.get(method, method),
                            'ont': ONT_LABEL.get(ont, ont),
                            'condition': COND_LABEL.get(cond, cond),
                            **row
                        })

    if not all_rows:
        print(f"No binned summaries found in {OUT_BASE}/")
        return

    fieldnames = ['method', 'ont', 'condition'] + [
        k for k in all_rows[0] if k not in ('method', 'ont', 'condition')]
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {out_path}")


# ── Print LaTeX flip-precision table to stdout ────────────────────────────────

def print_latex_table(results):
    print("\n% --- LaTeX flip-precision table (paste into thesis) ---")
    print(r"\begin{table}[!h]")
    print(r"\centering")
    print(r"\caption{Per-term flip precision of the taxon-consistency solver on the timeset split. "
          r"Flip precision = fraction of changed terms where the solver removed a false positive "
          r"(rather than a true positive). \label{tab:term-centric-flip}}")
    print(r"\small")
    print(r"\begin{tabular}{llcccc}")
    print(r"\toprule")
    print(r"Sub & Method & \multicolumn{2}{c}{Taxon-guided} & \multicolumn{2}{c}{Taxon-unknown} \\")
    print(r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}")
    print(r" & & Changed & Flip prec. & Changed & Flip prec. \\")
    print(r"\midrule")

    for ont in ONTS:
        first = True
        for method in METHODS:
            row_parts = [ONT_LABEL[ont] if first else "", METHOD_LABEL[method]]
            first = False
            for cond in CONDS:
                r = results.get((method, ont, cond), {})
                n = r.get('n_changed', '--')
                fp = r.get('flip_precision', '--')
                try:
                    fp_fmt = f"{float(fp):.2f}"
                except (TypeError, ValueError):
                    fp_fmt = '--'
                row_parts += [str(n), fp_fmt]
            print(" & ".join(row_parts) + r" \\")
        print(r"\midrule")

    print(r"\bottomrule")
    print(r"\end{tabular}")
    print(r"\end{table}")


if __name__ == "__main__":
    if not os.path.exists(LOG_FILE):
        sys.exit(f"Log file not found: {LOG_FILE}\nRun analysis/run_term_centric_timeset.sh first.")

    results = parse_log(LOG_FILE)
    print(f"Parsed {len(results)} runs from {LOG_FILE}")

    flip_path = os.path.join(OUT_BASE, "flip_precision_table.tsv")
    write_flip_precision_table(results, flip_path)

    binned_path = os.path.join(OUT_BASE, "binned_summary_combined.tsv")
    write_binned_combined(binned_path)

    print_latex_table(results)
