"""
Paired statistics on per-organism CAFA metrics before vs after taxon-guided
adjustment. For each metric, pairs organisms present in both the 'before'
(predictions) and 'after' (optimized) evaluation directories, then reports:
  - n paired organisms, median before/after, median delta
  - Wilcoxon signed-rank statistic W and two-sided p
  - rank-biserial effect size
  - organism-level bootstrap 95% CI on the mean delta

Eval JSONs are evaluate_directory.py output: evaluation_<fold>_<taxon>.json with
d['prediction_metrics'][subont][metric].
"""
import argparse
import glob
import json
import os
import re

import numpy as np
from scipy.stats import wilcoxon

ORG = re.compile(r"evaluation_(\d+)_(\d+)\.json$")


def load(dirpath, subont, metric):
    out = {}
    for f in glob.glob(os.path.join(dirpath, "evaluation_*_*.json")):
        m = ORG.search(os.path.basename(f))
        if not m:
            continue
        key = f"{m.group(1)}_{m.group(2)}"
        try:
            d = json.load(open(f))
            v = d.get("prediction_metrics", {}).get(subont, {}).get(metric)
            if v is not None:
                out[key] = float(v)
        except Exception:
            continue
    return out


def boot_ci(deltas, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(deltas), size=(n, len(deltas)))
    means = deltas[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before_dir", required=True)
    ap.add_argument("--after_dir", required=True)
    ap.add_argument("--subont", default="cc")
    ap.add_argument("--label", default="")
    ap.add_argument("--metrics", nargs="+", default=["fmax", "aupr", "avg_auc", "smin"])
    args = ap.parse_args()

    report = {"label": args.label, "subont": args.subont, "metrics": {}}
    for metric in args.metrics:
        b = load(args.before_dir, args.subont, metric)
        a = load(args.after_dir, args.subont, metric)
        keys = sorted(set(b) & set(a))
        if not keys:
            report["metrics"][metric] = {"n": 0}
            continue
        bv = np.array([b[k] for k in keys])
        av = np.array([a[k] for k in keys])
        d = av - bv
        nz = d[d != 0]
        if len(nz) >= 1:
            try:
                W, p = wilcoxon(av, bv, zero_method="wilcox", alternative="two-sided")
            except ValueError:
                W, p = float("nan"), 1.0
        else:
            W, p = float("nan"), 1.0
        npos, nneg = int((d > 0).sum()), int((d < 0).sum())
        rbc = (npos - nneg) / len(nz) if len(nz) else 0.0  # rank-biserial (sign-based)
        lo, hi = boot_ci(d) if len(d) > 1 else (float(d[0]), float(d[0]))
        report["metrics"][metric] = {
            "n": len(keys), "n_changed": int(len(nz)),
            "median_before": round(float(np.median(bv)), 4),
            "median_after": round(float(np.median(av)), 4),
            "mean_delta": round(float(d.mean()), 5),
            "median_delta": round(float(np.median(d)), 5),
            "n_improved": npos, "n_worsened": nneg,
            "wilcoxon_W": (round(float(W), 2) if W == W else None),
            "wilcoxon_p": (float(f"{p:.3e}") if p == p else None),
            "rank_biserial": round(float(rbc), 3),
            "boot95_mean_delta": [round(lo, 5), round(hi, 5)],
        }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
