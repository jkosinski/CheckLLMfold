"""Figures from results/summary.tsv (run with the data_analysis conda env, which has matplotlib).
fig_contact_maps.png: per target x replicate, reference contacts (lower triangle, grey) vs Fable's top-L ranked pairs
(upper triangle; blue correct, orange wrong), |i - j| >= 6.
fig_summary.png: top-L precision and TM-score per target, one dot per replicate, with chance and ceiling marks."""
import os, sys, json, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from structutil import parse, cb_ca, contact_map
from sscomp import dssp_lite
from evaluate import last_json

HIT, MISS, REF, INK, MUTED = '#2a78d6', '#eb6834', '#c9c8c2', '#0b0b0b', '#6e6d68'
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})

LABEL = {'ubiquitin_1UBQ': 'Ubiquitin, 1UBQ\n76 aa, in PDB',
         'cytb562_256B': 'Cytochrome b562, 256B\n106 aa all-α, in PDB',
         'A0A916UKZ0': 'A0A916UKZ0 (DUF2489)\n108 aa all-α, not in PDB',
         'fn3_tenascin_1TEN': 'Tenascin FN3, 1TEN\n89 aa all-β, in PDB',
         'A0A8J2YCU1': 'A0A8J2YCU1 (DUF2958)\n89 aa all-β, not in PDB',
         'thioredoxin_2TRX': 'Thioredoxin, 2TRX\n108 aa α/β, in PDB',
         'A0A5D0P002': 'A0A5D0P002 (DUF4288)\n107 aa α/β, not in PDB'}
SHORT = {'ubiquitin_1UBQ': 'Ubiquitin\n76 aa', 'cytb562_256B': 'Cyt\nb562', 'A0A916UKZ0': 'DUF\n2489',
         'fn3_tenascin_1TEN': 'Tenascin\nFN3', 'A0A8J2YCU1': 'DUF\n2958', 'thioredoxin_2TRX': 'Thio-\nredoxin',
         'A0A5D0P002': 'DUF\n4288'}
ORDER = ['ubiquitin_1UBQ', 'cytb562_256B', 'A0A916UKZ0', 'fn3_tenascin_1TEN', 'A0A8J2YCU1', 'thioredoxin_2TRX',
         'A0A5D0P002']
DARK = {'A0A916UKZ0', 'A0A8J2YCU1', 'A0A5D0P002'}
REFS = {'ubiquitin_1UBQ': ('data/1UBQ/1ubq.pdb', 'A'), 'cytb562_256B': ('data/256B/256b.pdb', 'A'),
        'fn3_tenascin_1TEN': ('data/1TEN/1ten.pdb', 'A'), 'thioredoxin_2TRX': ('data/2TRX/2trx.pdb', 'A')}
for a in DARK:
    REFS[a] = ('data/%s/%s.pdb' % (a, a), None)

rows = list(csv.DictReader(open(os.path.join(ROOT, 'results/summary.tsv')), delimiter='\t'))
rows = [r for r in rows if r['job'].startswith('job_')]


def ss_track(ax, ss, L, y, h):
    for k, c in enumerate(ss):
        if c == 'H':
            ax.add_patch(plt.Rectangle((k, y), 1, h, color='#e34948', lw=0))
        elif c == 'E':
            ax.add_patch(plt.Rectangle((k, y), 1, h, color='#6250d6', lw=0))


fig, axes = plt.subplots(len(ORDER), 3, figsize=(7.0, 2.25 * len(ORDER)))
for ri, t in enumerate(ORDER):
    ref, ch = REFS[t]
    res = [r for r in parse(os.path.join(ROOT, ref), ch) if 'CA' in r['atoms']]
    cb, ca = cb_ca(res)
    cm, _ = contact_map(cb)
    L = len(res)
    ss = dssp_lite(res)
    reps = sorted([r for r in rows if r['target'] == t], key=lambda r: r['rep'])
    for ci in range(3):
        ax = axes[ri, ci]
        ax.set_xlim(0, L); ax.set_ylim(L, -0.06 * L)
        ax.set_aspect('equal'); ax.set_xticks([1, L]); ax.set_yticks([1, L])
        if ci >= len(reps):
            ax.axis('off'); continue
        r = reps[ci]
        ans = last_json(open(os.path.join(ROOT, 'runs/jobs', r['job'], 'final.txt')).read())
        ranked, seen = [], set()
        for p in ans['pairs']:
            i, j = sorted((p['i'], p['j']))
            if j - i >= 6 and (i, j) not in seen and j <= L:
                seen.add((i, j)); ranked.append((i, j))
        ii, jj = np.where(np.triu(cm, 6))
        ax.scatter(ii + 0.5, jj + 0.5, s=2.2, marker='s', c=REF, lw=0)
        top = ranked[:L]
        hit = [(i, j) for i, j in top if cm[i - 1, j - 1]]
        miss = [(i, j) for i, j in top if not cm[i - 1, j - 1]]
        if miss:
            ax.scatter([j - 0.5 for i, j in miss], [i - 0.5 for i, j in miss], s=2.2, marker='s', c=MISS, lw=0)
        if hit:
            ax.scatter([j - 0.5 for i, j in hit], [i - 0.5 for i, j in hit], s=2.2, marker='s', c=HIT, lw=0)
        ax.plot([0, L], [0, L], color=MUTED, lw=0.4)
        ss_track(ax, ss, L, -0.05 * L, 0.03 * L)
        ax.set_title('rep %s: top-L P=%.2f, TM=%.2f' % (r['rep'], float(r['P_sep6_topL']), float(r['TM_best_hand'])),
                     fontsize=7, color=INK)
        if ci == 0:
            ax.set_ylabel(LABEL[t], fontsize=7.5)
fig.text(0.5, 0.005, 'Upper triangle: Fable 5.2 top-L ranked pairs (|i-j| ≥ 6), blue = in reference contact map, '
         'orange = not. Lower triangle: reference contacts (Cβ < 8 Å, |i-j| ≥ 6). Track: helix (red), strand (violet).',
         ha='center', fontsize=6.5, color=MUTED, wrap=True)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(os.path.join(ROOT, 'results/fig_contact_maps.png'), dpi=200)

INPDB, NOTPDB = '#6250d6', '#1baf7a'
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.4))
xs = []
for k, t in enumerate(ORDER):
    x = k + (k > 0) * 0.4 + ((k - 1) // 2) * 0.4 if k > 0 else 0
    xs.append(x)
    rr = [r for r in rows if r['target'] == t]
    if not rr:
        continue
    col = NOTPDB if t in DARK else INPDB
    p = [float(r['P_sep24_topL']) for r in rr]
    tm = [float(r['TM_best_hand']) for r in rr]
    jitter = np.linspace(-0.12, 0.12, len(rr)) if len(rr) > 1 else np.zeros(1)
    a1.scatter(x + jitter, p, s=24, c=col, zorder=3, edgecolor='white', lw=0.8)
    a2.scatter(x + jitter, tm, s=24, c=col, zorder=3, edgecolor='white', lw=0.8)
    a2.plot([x - 0.28, x + 0.28], [float(rr[0]['TM_ceiling_best_hand'])] * 2, color=INK, lw=1.4)
    # chance for long-range pairs: native LR contacts / all LR pairs
    L = int(rr[0]['L'])
    m = json.load(open(os.path.join(ROOT, 'runs/jobs', rr[0]['job'], 'metrics.json')))
    a1.plot([x - 0.28, x + 0.28], [m['native_contacts_sep24'] / ((L - 24) * (L - 23) / 2)] * 2, color=MUTED, lw=1.4)
for ax in (a1, a2):
    ax.set_xticks(xs); ax.set_xticklabels([SHORT[t] for t in ORDER], fontsize=6.5)
    ax.set_ylim(0, 1.06); ax.grid(axis='y', color='#e6e5e0', lw=0.6); ax.set_axisbelow(True)
    for (xa, xb), lab in (((xs[1], xs[2]), 'all-α, ~107 aa'), ((xs[3], xs[4]), 'all-β, 89 aa'), ((xs[5], xs[6]), 'α/β, ~107 aa')):
        ax.text((xa + xb) / 2, 1.03, lab, ha='center', fontsize=6.5, color=MUTED)
a1.set_title('Long-range precision, top-L pairs with |i-j| ≥ 24', fontsize=8, color=INK, loc='left')
a2.set_title('TM-score of model from final contacts (better hand)', fontsize=8, color=INK, loc='left')
from matplotlib.lines import Line2D
h = [Line2D([], [], marker='o', ls='', color=INPDB, label='in PDB, well studied'),
     Line2D([], [], marker='o', ls='', color=NOTPDB, label='no relative in PDB (AFDB model)'),
     Line2D([], [], color=MUTED, lw=1.4, label='chance (left)'),
     Line2D([], [], color=INK, lw=1.4, label='fold() on true contacts (right)')]
fig.legend(handles=h, loc='lower center', ncol=4, fontsize=6.5, frameon=False)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(os.path.join(ROOT, 'results/fig_summary.png'), dpi=200)
print('wrote figures')
