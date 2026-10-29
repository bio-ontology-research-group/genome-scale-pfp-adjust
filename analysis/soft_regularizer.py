#!/usr/bin/env python3
"""Historical greedy margin-threshold diagnostic.

This exploratory script is not an implementation of OCELOT, Semantic-Based
Regularization, or the revised jointly optimized soft objective. Its old
numerical conclusions are not accepted revision evidence. Use the controlled
soft-constraint conditions in analysis/revision_audit_run.py instead.
"""

import argparse
import glob
import os
import re

FN = re.compile(r"_fold_(\d+)_taxon_(\d+)\.tsv$")


def parse_pred(path):
    """protein -> {term: score} for entries with a score field."""
    out = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            d = {}
            for tok in parts[1:]:
                if "|" in tok:
                    t, s = tok.split("|", 1)
                    try:
                        d[t] = float(s)
                    except ValueError:
                        pass
            out[parts[0]] = d
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--predictions_dir", required=True)
    ap.add_argument("--adjusted_dir", required=True,
                    help="With-taxon Stage-1 output dir (optimized_fold_*_taxon_*.tsv).")
    ap.add_argument("--threshold", type=float, default=0.3)
    ap.add_argument("--label", default="run")
    ap.add_argument("--lambdas", default="0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    tau = args.threshold
    lambdas = [float(x) for x in args.lambdas.split(",")]

    margins = []          # margin above tau for every hard-removed (violating) term
    n_org = 0
    org_has_violation = 0
    for adj in sorted(glob.glob(os.path.join(args.adjusted_dir, "optimized_fold_*_taxon_*.tsv"))):
        m = FN.search(os.path.basename(adj))
        if not m:
            continue
        pred = os.path.join(args.predictions_dir,
                            f"predictions_fold_{m.group(1)}_taxon_{m.group(2)}.tsv")
        if not os.path.exists(pred):
            continue
        n_org += 1
        P, A = parse_pred(pred), parse_pred(adj)
        org_viol = 0
        for prot, terms in P.items():
            adj_terms = A.get(prot, {})
            for t, s in terms.items():
                if s > tau:  # predicted above threshold
                    a = adj_terms.get(t, 0.0)
                    if a <= tau:  # hard solver removed it -> taxon violation
                        margins.append(s - tau)
                        org_viol += 1
        if org_viol:
            org_has_violation += 1

    total = len(margins)
    rows = []
    # Hard solver: removes all violations -> 0 residual, full consistency.
    rows.append(("hard(exact)", 0, 100.0, sum(margins)))
    for lam in lambdas:
        # Soft penalty lambda removes a violation only if its margin < lambda.
        removed = [mg for mg in margins if mg < lam]
        residual = total - len(removed)
        pct = 100.0 * len(removed) / total if total else 100.0
        rows.append((f"soft(lambda={lam:g})", residual, pct, sum(removed)))

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(f"# soft-regularizer vs hard taxon solver | {args.label} | "
                 f"tau={tau} | organisms={n_org} | with-violation={org_has_violation} | "
                 f"total_violating_terms={total}\n")
        fh.write("mode\tresidual_violations\tpct_removed\tremoval_cost\n")
        for name, residual, pct, cost in rows:
            fh.write(f"{name}\t{residual}\t{pct:.1f}\t{cost:.3f}\n")

    print(f"organisms={n_org} with-violation={org_has_violation} total_violating_terms={total}")
    for name, residual, pct, cost in rows:
        print(f"  {name:18s} residual={residual:5d} removed={pct:5.1f}% cost={cost:.3f}")
    print(f"wrote {args.output}")


if __name__ == '__main__':
    raise SystemExit('This historical evaluator is retired. Use the documented reproducibility/revision workflow.')
