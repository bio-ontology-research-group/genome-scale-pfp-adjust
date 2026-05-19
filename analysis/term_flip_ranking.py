"""
Term-centric flip ranking.

For each (GO term) across the timeset test proteins, count the number of flips
the solver applies and the total "cost" of those flips (sum of absolute score
differences), restricted to flips whose *original* prediction score exceeds a
threshold (default 0.1). Flips below the threshold are considered too weak to
be biologically interesting.

A "flip" on a (protein, term) pair is any non-zero difference between the
original prediction score and the post-solver score. Because the solver is
removal-only for taxon consistency, flips on that stage are all removals;
complex coherence can in principle add terms as well, so we keep the
direction (removed / added) as an output column for generality.

Outputs:
  - <out>_by_count.tsv  : top-K terms ranked by flip count
  - <out>_by_cost.tsv   : top-K terms ranked by cumulative flip cost
  - <out>_full.tsv      : per-term full table (every term with >=1 flip)
  - <out>_summary.tsv   : global totals (n_terms, n_flips, n_proteins ...)
"""

from typing import Dict, List, Optional, Tuple
from collections import defaultdict
import argparse
import glob
import os
import sys


def load_test_proteins(test_proteins_file: str) -> Dict[str, List[str]]:
    test_proteins: Dict[str, List[str]] = defaultdict(list)
    with open(test_proteins_file, 'r') as f:
        header = f.readline().strip().split('\t')
        taxon_id_index = header.index('taxon_id')
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            protein_id = parts[0]
            if protein_id == 'protein_id':
                continue
            taxon_id = parts[taxon_id_index]
            if taxon_id == 'NA':
                continue
            test_proteins[taxon_id].append(protein_id)
    return test_proteins


def load_per_taxon_tsv(directory: str, pattern: str,
                       test_proteins: Dict[str, List[str]]
                       ) -> Dict[str, Dict[str, float]]:
    """Load {protein_id -> {go_id -> score}} from per-taxon TSV files.

    `pattern` is a glob pattern accepting one `{taxon_id}` placeholder.
    """
    out: Dict[str, Dict[str, float]] = {}
    for taxon_id, protein_ids in test_proteins.items():
        glob_pat = pattern.format(taxon_id=taxon_id)
        files = glob.glob(os.path.join(directory, glob_pat))
        if not files:
            continue
        protein_id_set = set(protein_ids)
        with open(files[0], 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if not parts:
                    continue
                pid = parts[0]
                if pid not in protein_id_set:
                    continue
                term_scores: Dict[str, float] = {}
                for pred in parts[1:]:
                    if '|' not in pred:
                        continue
                    go_id, score_str = pred.split('|', 1)
                    try:
                        score = float(score_str)
                    except ValueError:
                        continue
                    term_scores[go_id] = score
                out[pid] = term_scores
    return out


def load_go_metadata(go_file: str, subontology: str
                     ) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Return (go_id -> name, go_id -> depth) for terms in the requested subontology.

    Depth = shortest path length from the subontology root; root has depth 0.
    Terms not in the requested subontology are omitted.
    """
    names: Dict[str, str] = {}
    parents: Dict[str, List[str]] = defaultdict(list)
    namespace: Dict[str, str] = {}
    obsolete: Dict[str, bool] = {}

    ns_map = {
        'cc': 'cellular_component',
        'mf': 'molecular_function',
        'bp': 'biological_process',
    }
    target_ns = ns_map.get(subontology.lower())

    current_id: Optional[str] = None
    in_term = False
    with open(go_file, 'r') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if line == '[Term]':
                in_term = True
                current_id = None
                continue
            if line.startswith('['):
                in_term = False
                current_id = None
                continue
            if not in_term:
                continue
            if line.startswith('id: '):
                current_id = line[4:].strip()
            elif line.startswith('name: ') and current_id:
                names[current_id] = line[6:].strip()
            elif line.startswith('namespace: ') and current_id:
                namespace[current_id] = line[11:].strip()
            elif line.startswith('is_obsolete: true') and current_id:
                obsolete[current_id] = True
            elif line.startswith('is_a: ') and current_id:
                parent = line[6:].split(' !', 1)[0].strip()
                parents[current_id].append(parent)
            elif line.startswith('relationship: part_of ') and current_id:
                parent = line[len('relationship: part_of '):].split(' !', 1)[0].strip()
                parents[current_id].append(parent)

    if target_ns is not None:
        kept = {gid for gid, ns in namespace.items() if ns == target_ns and not obsolete.get(gid, False)}
    else:
        kept = set(names.keys()) - {gid for gid, obs in obsolete.items() if obs}

    depth: Dict[str, int] = {}
    def compute_depth(gid: str, seen: set) -> int:
        if gid in depth:
            return depth[gid]
        if gid in seen:
            return 0
        ps = [p for p in parents.get(gid, []) if p in kept]
        if not ps:
            depth[gid] = 0
            return 0
        d = 1 + min(compute_depth(p, seen | {gid}) for p in ps)
        depth[gid] = d
        return d

    for gid in kept:
        compute_depth(gid, set())

    names_out = {gid: names.get(gid, '') for gid in kept}
    return names_out, depth


def load_taxon_constraints(constraints_file: str) -> Dict[str, str]:
    """Return go_id -> 'only_in' / 'never_in' / 'both'. Missing = no constraint."""
    out: Dict[str, str] = {}
    if not constraints_file or not os.path.exists(constraints_file):
        return out
    with open(constraints_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            gid, rel = parts[0], parts[1]
            rel_key = 'only_in' if 'only' in rel.lower() else ('never_in' if 'never' in rel.lower() else rel)
            if gid in out and out[gid] != rel_key:
                out[gid] = 'both'
            else:
                out[gid] = rel_key
    return out


def load_ground_truth(annotations_dir: str, subontology: str,
                      test_proteins: Dict[str, List[str]]
                      ) -> Dict[str, set]:
    """Optional: load ground-truth (protein -> set of GO terms) for flip-precision.

    Expects per-taxon annotation TSVs keyed by taxon_id; schema may vary between
    projects, so we only use this if the directory exists.
    Returns {} if annotations_dir is missing.
    """
    if not annotations_dir or not os.path.isdir(annotations_dir):
        return {}
    gt: Dict[str, set] = {}
    ns_map = {
        'cc': 'C', 'mf': 'F', 'bp': 'P',
    }
    target_aspect = ns_map.get(subontology.lower())
    for taxon_id, protein_ids in test_proteins.items():
        files = glob.glob(os.path.join(annotations_dir, f'*_taxon_{taxon_id}.*'))
        if not files:
            files = glob.glob(os.path.join(annotations_dir, f'{taxon_id}.*'))
        if not files:
            continue
        protein_id_set = set(protein_ids)
        for fp in files:
            try:
                with open(fp, 'r') as f:
                    for line in f:
                        if line.startswith('!'):
                            continue
                        parts = line.rstrip('\n').split('\t')
                        if len(parts) < 9:
                            continue
                        pid = parts[1] if parts[1] else parts[0]
                        if pid not in protein_id_set:
                            continue
                        go_id = parts[4]
                        aspect = parts[8] if len(parts) > 8 else ''
                        if target_aspect and aspect != target_aspect:
                            continue
                        gt.setdefault(pid, set()).add(go_id)
            except Exception:
                continue
    return gt


def rank_flips(
    pre: Dict[str, Dict[str, float]],
    post: Dict[str, Dict[str, float]],
    score_threshold: float,
) -> Tuple[Dict[str, dict], dict]:
    """Compute per-term flip stats.

    A flip = score changed between pre and post.
    Only flips where pre-score > threshold (for removals) OR post-score > threshold
    (for additions) are counted.

    Returns (per_term_stats, global_stats).
    per_term_stats[go_id] = {
        n_flips, n_removals, n_additions,
        cost, cost_removals, cost_additions,
        n_proteins, mean_pre_score, max_pre_score,
    }
    """
    proteins = set(pre.keys()) | set(post.keys())
    stats: Dict[str, dict] = defaultdict(lambda: {
        'n_flips': 0, 'n_removals': 0, 'n_additions': 0,
        'cost': 0.0, 'cost_removals': 0.0, 'cost_additions': 0.0,
        'proteins': set(),
        'pre_scores': [],
    })

    total_flips = 0
    total_cost = 0.0
    total_pre_terms_above = 0
    for pid in proteins:
        pre_ts = pre.get(pid, {})
        post_ts = post.get(pid, {})
        for go_id, pre_score in pre_ts.items():
            if pre_score > score_threshold:
                total_pre_terms_above += 1
            post_score = post_ts.get(go_id, 0.0)
            if pre_score == post_score:
                continue
            if pre_score > score_threshold and post_score < pre_score:
                delta = pre_score - post_score
                s = stats[go_id]
                s['n_flips'] += 1
                s['n_removals'] += 1
                s['cost'] += delta
                s['cost_removals'] += delta
                s['proteins'].add(pid)
                s['pre_scores'].append(pre_score)
                total_flips += 1
                total_cost += delta
        for go_id, post_score in post_ts.items():
            if go_id in pre_ts:
                continue
            if post_score > score_threshold:
                s = stats[go_id]
                s['n_flips'] += 1
                s['n_additions'] += 1
                s['cost'] += post_score
                s['cost_additions'] += post_score
                s['proteins'].add(pid)
                s['pre_scores'].append(0.0)
                total_flips += 1
                total_cost += post_score

    per_term: Dict[str, dict] = {}
    for go_id, s in stats.items():
        scores = s['pre_scores']
        per_term[go_id] = {
            'n_flips': s['n_flips'],
            'n_removals': s['n_removals'],
            'n_additions': s['n_additions'],
            'cost': s['cost'],
            'cost_removals': s['cost_removals'],
            'cost_additions': s['cost_additions'],
            'n_proteins': len(s['proteins']),
            'mean_pre_score': sum(scores) / len(scores) if scores else 0.0,
            'max_pre_score': max(scores) if scores else 0.0,
        }

    global_stats = {
        'n_proteins_evaluated': len(proteins),
        'n_terms_with_flips': len(per_term),
        'n_flips_total': total_flips,
        'cost_total': total_cost,
        'n_pre_terms_above_threshold': total_pre_terms_above,
    }
    return per_term, global_stats


def write_table(rows: List[dict], columns: List[str], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    with open(path, 'w') as f:
        f.write('\t'.join(columns) + '\n')
        for r in rows:
            f.write('\t'.join(
                (f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, '')))
                for c in columns
            ) + '\n')


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--predictions_dir', required=True,
                   help='Directory with predictions_fold_*_taxon_*.tsv')
    p.add_argument('--optimized_dir', required=True,
                   help='Directory with optimized_*_taxon_*.tsv')
    p.add_argument('--test_proteins_file', required=True,
                   help='TSV with protein_id + taxon_id columns for test set')
    p.add_argument('--go_file', required=True, help='go-basic.obo')
    p.add_argument('--subontology', required=True, choices=['cc', 'mf', 'bp'])
    p.add_argument('--output_prefix', required=True,
                   help='Prefix for output files; _by_count.tsv, _by_cost.tsv, '
                        '_full.tsv, _summary.tsv are appended.')
    p.add_argument('--score_threshold', type=float, default=0.1,
                   help='Minimum ORIGINAL prediction score for a flip to be counted '
                        '(default: 0.1).')
    p.add_argument('--top_k', type=int, default=25)
    p.add_argument('--taxon_constraints_file', default=None,
                   help='Optional TSV: go_id \\t only_in_taxon|never_in_taxon \\t taxon_id')
    p.add_argument('--annotations_dir', default=None,
                   help='Optional ground-truth annotations dir for flip-precision.')
    args = p.parse_args(argv)

    print(f"[load] test proteins: {args.test_proteins_file}", file=sys.stderr)
    test_proteins = load_test_proteins(args.test_proteins_file)
    n_proteins = sum(len(v) for v in test_proteins.values())
    print(f"[load]   -> {len(test_proteins)} taxa, {n_proteins} proteins", file=sys.stderr)

    print(f"[load] predictions: {args.predictions_dir}", file=sys.stderr)
    pre = load_per_taxon_tsv(args.predictions_dir,
                             'predictions_fold_*_taxon_{taxon_id}.tsv',
                             test_proteins)
    print(f"[load]   -> {len(pre)} proteins with predictions", file=sys.stderr)

    print(f"[load] optimized: {args.optimized_dir}", file=sys.stderr)
    post = load_per_taxon_tsv(args.optimized_dir,
                              'optimized_*_taxon_{taxon_id}.tsv',
                              test_proteins)
    print(f"[load]   -> {len(post)} proteins with optimized output", file=sys.stderr)

    print(f"[load] GO metadata: {args.go_file}", file=sys.stderr)
    go_names, go_depth = load_go_metadata(args.go_file, args.subontology)

    constraints_map = load_taxon_constraints(args.taxon_constraints_file) if args.taxon_constraints_file else {}
    gt = load_ground_truth(args.annotations_dir, args.subontology, test_proteins) if args.annotations_dir else {}

    print(f"[rank] threshold={args.score_threshold}", file=sys.stderr)
    per_term, global_stats = rank_flips(pre, post, args.score_threshold)

    rows: List[dict] = []
    for go_id, s in per_term.items():
        row = {
            'go_id': go_id,
            'name': go_names.get(go_id, ''),
            'depth': go_depth.get(go_id, -1),
            'taxon_constraint': constraints_map.get(go_id, ''),
            **s,
        }
        if gt:
            tp = 0; fp = 0
            # Compute flip-precision: fraction of flipped (protein, term) pairs
            # that were NOT in ground truth (i.e. removing a FP is "correct").
            # This requires tracking which proteins flipped for each term; we
            # only have counts here, so we leave flip_precision blank unless we
            # re-walk the data. (Kept as placeholder for future extension.)
            row['flip_precision'] = ''
        rows.append(row)

    top_columns = ['go_id', 'name', 'depth', 'taxon_constraint',
                   'n_flips', 'n_removals', 'n_additions',
                   'cost', 'cost_removals', 'cost_additions',
                   'n_proteins', 'mean_pre_score', 'max_pre_score']

    by_count = sorted(rows, key=lambda r: (-r['n_flips'], -r['cost']))[:args.top_k]
    by_cost = sorted(rows, key=lambda r: (-r['cost'], -r['n_flips']))[:args.top_k]
    full = sorted(rows, key=lambda r: (-r['n_flips'], -r['cost']))

    write_table(by_count, top_columns, f'{args.output_prefix}_by_count.tsv')
    write_table(by_cost, top_columns, f'{args.output_prefix}_by_cost.tsv')
    write_table(full, top_columns, f'{args.output_prefix}_full.tsv')

    with open(f'{args.output_prefix}_summary.tsv', 'w') as f:
        f.write('key\tvalue\n')
        for k, v in global_stats.items():
            f.write(f'{k}\t{v}\n')
        f.write(f'score_threshold\t{args.score_threshold}\n')
        f.write(f'subontology\t{args.subontology}\n')

    print(f"[done] terms with flips: {global_stats['n_terms_with_flips']}, "
          f"total flips: {global_stats['n_flips_total']}, "
          f"total cost: {global_stats['cost_total']:.2f}",
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
