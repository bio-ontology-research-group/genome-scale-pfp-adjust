"""
Compute the Stage-1 (taxon) flip rate between a predictions directory and an
adjusted (optimized) directory at a given threshold tau, aggregated across all
matching per-organism files. Stage 1 is demotion-only, so a "flip" is a
(protein, GO term) annotation that is above tau in the prediction but absent
(or below tau) in the adjusted output.

Outputs one summary line:  preds_above_tau  flips  flip_rate_pct
"""
import argparse
import glob
import os
import re

FN = re.compile(r"_fold_(\d+)_taxon_(\d+)\.tsv$")


def above(path, tau):
    """protein -> set(terms with score > tau)."""
    d = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if not p or not p[0]:
                continue
            s = set()
            for tok in p[1:]:
                if "|" in tok:
                    go, sc = tok.split("|", 1)
                    try:
                        if float(sc) > tau:
                            s.add(go)
                    except ValueError:
                        pass
            d[p[0]] = s
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions_dir", required=True)
    ap.add_argument("--optimized_dir", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    preds_n = flips = 0
    npairs = 0
    for opt in sorted(glob.glob(os.path.join(args.optimized_dir, "optimized_fold_*_taxon_*.tsv"))):
        m = FN.search(os.path.basename(opt))
        if not m:
            continue
        fold, taxon = m.group(1), m.group(2)
        pred = os.path.join(args.predictions_dir, f"predictions_fold_{fold}_taxon_{taxon}.tsv")
        if not os.path.exists(pred):
            continue
        npairs += 1
        tp = above(pred, args.threshold)
        to = above(opt, args.threshold)
        for prot, pset in tp.items():
            preds_n += len(pset)
            flips += len(pset - to.get(prot, set()))
    rate = (100.0 * flips / preds_n) if preds_n else 0.0
    print(f"{args.label}\ttau={args.threshold}\torganisms={npairs}\tpreds_above_tau={preds_n}\tflips={flips}\tflip_rate_pct={rate:.4f}")


if __name__ == "__main__":
    main()
