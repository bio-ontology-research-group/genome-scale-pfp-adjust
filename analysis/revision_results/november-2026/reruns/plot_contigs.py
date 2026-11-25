"""Stacked-bar figures for the contig lineage rerun (Supplementary Figures S1, S2)."""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RES = Path(sys.argv[1])
OUT = Path(sys.argv[2])
ORGS = [('83332', 'M. tuberculosis H37Rv'), ('3572', 'Saponaria officinalis'), ('69014', 'T. kodakarensis')]
PRED = {'seq_sim_timeset': 'Seq-Sim', 'deepgo-se': 'DeepGO-SE', 'sprof-go': 'SPROF-GO'}
LENGTHS = ['1000', '10000', '100000', '1000000']
LLAB = ['1 kb', '10 kb', '100 kb', '1 Mb']
CATS = [('no_annotations', 'no annotations above $\\tau$', '#d9d9d9'),
        ('unresolved', 'no required taxon', '#9e9e9e'),
        ('wrong', 'wrong lineage', '#c0504d'),
        ('correct_domain', 'correct: domain', '#c6dbef'),
        ('correct_phylum', 'correct: phylum', '#6baed6'),
        ('correct_deeper', 'correct: class or deeper', '#08519c'),
        ('other', 'not certified / error', '#fdae6b')]


def cat(c):
    if c in {'correct_class', 'correct_order', 'correct_family', 'correct_genus', 'correct_species'}:
        return 'correct_deeper'
    if c in {'not_certified', 'solver_error', 'correct_none'}:
        return 'other'
    return c


def panel(ax, data, title):
    bottom = [0] * len(LENGTHS)
    for key, label, colour in CATS:
        vals = []
        for L in LENGTHS:
            s = data['summary'].get(L, {})
            tot = sum(s.values()) or 1
            vals.append(100 * sum(v for k, v in s.items() if cat(k) == key) / tot)
        ax.bar(range(len(LENGTHS)), vals, bottom=bottom, color=colour, label=label, width=.75)
        bottom = [b + v for b, v in zip(bottom, vals)]
    ax.set_xticks(range(len(LENGTHS)), LLAB, fontsize=8)
    ax.set_ylim(0, 100)
    ax.set_title(title, fontsize=9)


def figure(cols, name):
    fig, axes = plt.subplots(3, len(cols), figsize=(2.6 * len(cols), 7), sharey=True)
    for r, (tax, oname) in enumerate(ORGS):
        for c, (pred, asp, ctitle) in enumerate(cols):
            data = json.loads((RES / f'{tax}_{pred}_{asp}.json').read_text())
            panel(axes[r][c], data, ctitle if r == 0 else '')
            if c == 0:
                axes[r][c].set_ylabel(f'{oname}\n% of contigs', fontsize=9, style='italic')
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, .07, 1, 1))
    fig.savefig(OUT / name)


figure([(p, 'bp', n) for p, n in PRED.items()], 'metagenomics_fig_main_condense.pdf')
figure([('seq_sim_timeset', a, a.upper()) for a in ['cc', 'mf', 'bp']], 'metagenomics_fig_summary.pdf')
