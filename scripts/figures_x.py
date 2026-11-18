"""Standalone 16:9 figures (1600x900) for an X thread. Run with the data_analysis conda env (matplotlib).
Writes results/x_thread/x1_control_ubiquitin.png ... x4_summary.png"""
import os, sys, csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from structutil import parse, cb_ca, contact_map
from evaluate import last_json

OUT = os.path.join(ROOT, 'results/x_thread')
os.makedirs(OUT, exist_ok=True)
HIT, MISS, REF, INK, MUTED, GRID, BG = '#2a78d6', '#eb6834', '#cfcec8', '#0b0b0b', '#5f5e5a', '#e6e5e0', '#ffffff'
INPDB, NOTPDB = '#6250d6', '#1baf7a'
plt.rcParams.update({'font.size': 11, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': BG, 'savefig.facecolor': BG, 'font.family': 'DejaVu Sans'})
FOOT = ('Claude Fable 5.2 (effort high), Sergey Ovchinnikov\'s contact-map prompt, run headless with no network.\n'
        'Contact: Cβ–Cβ < 8 Å. L = protein length.')

rows = list(csv.DictReader(open(os.path.join(ROOT, 'results/summary.tsv')), delimiter='\t'))
REFS = {'ubiquitin_1UBQ': ('data/1UBQ/1ubq.pdb', 'A'), 'cytb562_256B': ('data/256B/256b.pdb', 'A'),
        'fn3_tenascin_1TEN': ('data/1TEN/1ten.pdb', 'A'), 'thioredoxin_2TRX': ('data/2TRX/2trx.pdb', 'A')}
for a in ('A0A916UKZ0', 'A0A8J2YCU1', 'A0A5D0P002'):
    REFS[a] = ('data/%s/%s.pdb' % (a, a), None)


def refmap(t):
    f, ch = REFS[t]
    res = [r for r in parse(os.path.join(ROOT, f), ch) if 'CA' in r['atoms']]
    cb, _ = cb_ca(res)
    return contact_map(cb)[0], len(res)


def ranked(run_dir, L):
    ans = last_json(open(os.path.join(run_dir, 'final.txt')).read())
    out, seen = [], set()
    for p in ans['pairs']:
        i, j = sorted((p['i'], p['j']))
        if j - i >= 6 and j <= L and (i, j) not in seen:
            seen.add((i, j)); out.append((i, j))
    return out


def cmap_panel(ax, cm, L, pairs, title, sub=None, title_color=INK):
    ms = max(1.2, 9000 / L ** 2)
    ii, jj = np.where(np.triu(cm, 6))
    ax.scatter(ii + 0.5, jj + 0.5, s=ms, marker='s', c=REF, lw=0)
    top = pairs[:L]
    for sel, col in (([p for p in top if not cm[p[0] - 1, p[1] - 1]], MISS), ([p for p in top if cm[p[0] - 1, p[1] - 1]], HIT)):
        if sel:
            ax.scatter([j - 0.5 for i, j in sel], [i - 0.5 for i, j in sel], s=ms, marker='s', c=col, lw=0)
    ax.plot([0, L], [0, L], color=MUTED, lw=0.5)
    ax.set_xlim(0, L); ax.set_ylim(L, 0); ax.set_aspect('equal')
    ax.set_xticks([1, L]); ax.set_yticks([1, L]); ax.tick_params(labelsize=8, length=2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(True)
    ax.set_title(title, fontsize=11, color=title_color, pad=4)
    if sub:
        ax.text(0.5, -0.13, sub, transform=ax.transAxes, ha='center', va='top', fontsize=9.5, color=INK)


def legend_maps(fig, y):
    h = [Patch(color=HIT, label='Fable top-L pair, real contact'), Patch(color=MISS, label='Fable top-L pair, not a contact'),
         Patch(color=REF, label='contacts in the reference structure')]
    fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, y), ncol=3, fontsize=9.5, frameon=False,
               handlelength=1.0, handleheight=1.0)


def frac(cm, pairs, L, minsep=6):
    sel = [p for p in pairs if p[1] - p[0] >= minsep][:L]
    return np.mean([cm[i - 1, j - 1] for i, j in sel]) if sel else float('nan')


def foot(fig, text=FOOT):
    fig.text(0.012, 0.012, text, fontsize=7.5, color=MUTED, ha='left', va='bottom')


# ---------- 1. ubiquitin control ----------
cm, L = refmap('ubiquitin_1UBQ')
fig = plt.figure(figsize=(8, 4.5), dpi=200)
fig.suptitle('Ubiquitin (PDB 1UBQ, 76 aa): Fable 5.2 contact predictions', fontsize=14, color=INK, x=0.03, ha='left', y=0.965)
fig.text(0.03, 0.865, 'Upper triangle: Fable\'s 76 top-ranked pairs (|i−j| ≥ 6). Lower triangle: crystal structure.',
         fontsize=9.5, color=MUTED, ha='left')
runs = [r for r in rows if r['target'] == 'ubiquitin_1UBQ']
runs = sorted([r for r in runs if r['job'].startswith('job_')], key=lambda r: r['rep']) + [r for r in runs if not r['job'].startswith('job_')]
runs = [r for r in runs if os.path.exists(os.path.join(ROOT, 'runs/jobs' if r['job'].startswith('job_') else 'runs', r['job'], 'final.txt'))]
for k, r in enumerate(runs):
    ax = fig.add_axes([0.04 + k * 0.245, 0.33, 0.205, 0.205 * 16 / 9])
    d = os.path.join(ROOT, 'runs/jobs', r['job']) if r['job'].startswith('job_') else os.path.join(ROOT, 'runs', r['job'])
    pr = ranked(d, L)
    name = 'our run %s' % r['rep'] if r['job'].startswith('job_') else 'Sergey\'s chat'
    cmap_panel(ax, cm, L, pr, name, 'top-L precision %.2f\nlong-range (|i−j| ≥ 24) %.2f' % (frac(cm, pr, L), frac(cm, pr, L, 24)))
legend_maps(fig, 0.075)
foot(fig, 'Claude Fable 5.2, Sergey Ovchinnikov\'s prompt. Our runs: headless, effort high, no network.\n'
     'Chance precision 0.06. Contact: Cβ–Cβ < 8 Å.')
fig.savefig(os.path.join(OUT, 'x1_control_ubiquitin.png'))
plt.close(fig)

# ---------- 2. selection of proteins with no PDB relative ----------
fig = plt.figure(figsize=(8, 4.5), dpi=200)
fig.suptitle('Choosing proteins with no relative in the PDB', fontsize=14, color=INK, x=0.03, ha='left', y=0.965)
ax = fig.add_axes([0.03, 0.14, 0.44, 0.70]); ax.axis('off')
steps = [(500, 'bacterial DUF proteins, 75–110 aa,\nAlphaFold DB model, not in PDB'),
         (79, 'mean pLDDT ≥ 90'),
         (78, 'no PDB sequence hit\n(RCSB MMseqs2, E < 1)'),
         (27, 'compact; Foldseek search\nagainst PDB100'),
         (6, 'no Foldseek hit, E < 0.5'),
         (3, 'Pfam family: no PDB\nstructure, no clan')]
mx = np.log10(500)
for k, (n, lab) in enumerate(steps):
    y = 1 - k / len(steps)
    w = 0.20 * np.log10(n) / mx + 0.03
    ax.add_patch(plt.Rectangle((0, y - 0.13), w, 0.11, color=NOTPDB if k == len(steps) - 1 else '#b9d7cb', lw=0))
    ax.text(w + 0.02, y - 0.075, str(n), fontsize=12, color=INK, va='center', fontweight='bold')
    ax.text(0.40, y - 0.075, lab, fontsize=8.6, color=INK, va='center')
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
tab = fig.add_axes([0.505, 0.14, 0.475, 0.70]); tab.axis('off')
cols = ['', 'DUF2489', 'DUF2958', 'DUF4288']
data = [['UniProt', 'A0A916UKZ0', 'A0A8J2YCU1', 'A0A5D0P002'],
        ['class, aa', 'all-α, 108', 'all-β, 89', 'α/β, 107'],
        ['mean pLDDT', '97.5', '91.1', '94.4'],
        ['pLDDT < 70', 'none', '3 at ends', '4 at ends'],
        ['mean PAE', '2.4 Å', '3.8 Å', '3.6 Å'],
        ['Foldseek E', '1.1', '0.52', '0.59'],
        ['best TM, PDB', '0.60', '0.52', '0.54']]
t = tab.table(cellText=data, colLabels=cols, loc='upper left', cellLoc='center', colLoc='center',
              colWidths=[0.25, 0.25, 0.25, 0.25])
t.auto_set_font_size(False); t.set_fontsize(8.6); t.scale(1, 1.75)
for (r, c), cell in t.get_celld().items():
    cell.set_edgecolor(GRID); cell.set_linewidth(0.6)
    if r == 0:
        cell.set_text_props(color=INK, fontweight='bold')
    if c == 0:
        cell.set_text_props(color=MUTED)
    if r == 1 and c > 0:
        cell.set_text_props(fontsize=7.6)
tab.text(0.0, 0.0, 'Best TM-score hits are unrelated proteins of similar shape\n(other helical bundles, an Sm-like barrel, a ribosomal protein).',
         fontsize=8, color=MUTED, transform=tab.transAxes, va='bottom')
foot(fig, 'Reference structures for these three: AlphaFold DB v6 models.\n'
     'Foldseek server, 3Di+AA and TM-align modes; TM-score normalised by query length.')
fig.savefig(os.path.join(OUT, 'x2_selection.png'))
plt.close(fig)

# ---------- 3. class-matched contact maps, one image per class ----------
pairs_cls = [('a', 'all-α, about 107 aa', 'cytb562_256B', 'Cytochrome b562 (PDB 256B)', 'A0A916UKZ0', 'DUF2489 (UniProt A0A916UKZ0)'),
             ('b', 'all-β, 89 aa', 'fn3_tenascin_1TEN', 'Tenascin FN3 domain (PDB 1TEN)', 'A0A8J2YCU1', 'DUF2958 (UniProt A0A8J2YCU1)'),
             ('c', 'α/β, about 107 aa', 'thioredoxin_2TRX', 'Thioredoxin (PDB 2TRX)', 'A0A5D0P002', 'DUF4288 (UniProt A0A5D0P002)')]
for tag, cls, famous, fname, dark, dname in pairs_cls:
    fig = plt.figure(figsize=(8, 4.5), dpi=200)
    fig.suptitle('Same size and class (%s): well studied vs no PDB relative' % cls, fontsize=13.5, color=INK, x=0.03, ha='left', y=0.965)
    fig.text(0.03, 0.875, 'Upper triangle: Fable\'s L top-ranked pairs (|i−j| ≥ 6). Lower triangle: reference structure.',
             fontsize=9.5, color=MUTED, ha='left')
    for k, (t, nm, col, lab) in enumerate(((famous, fname, INPDB, 'in PDB, well studied'), (dark, dname, NOTPDB, 'no relative in PDB'))):
        cm, L = refmap(t)
        cand = [r for r in rows if r['target'] == t and r['job'].startswith('job_')]
        best = max(cand, key=lambda r: float(r['P_sep24_topL']))
        pr = ranked(os.path.join(ROOT, 'runs/jobs', best['job']), L)
        ax = fig.add_axes([0.05 + k * 0.48, 0.215, 0.30, 0.5])
        cmap_panel(ax, cm, L, pr, '')
        ax.text(0.0, 1.17, lab, transform=ax.transAxes, fontsize=10, color=col, fontweight='bold', ha='left')
        ax.text(0.0, 1.05, nm, transform=ax.transAxes, fontsize=10.5, color=INK, ha='left')
        allp = [float(r['P_sep24_topL']) for r in cand]
        ax.text(1.05, 0.62, 'long-range\nprecision', transform=ax.transAxes, fontsize=9.5, color=MUTED, ha='left', va='bottom')
        ax.text(1.05, 0.6, '%.2f' % frac(cm, pr, L, 24), transform=ax.transAxes, fontsize=22, color=INK, ha='left', va='top')
        if len(cand) > 1:
            ax.text(1.05, 0.38, 'best of %d runs\nall: %s' % (len(cand), ', '.join('%.2f' % v for v in sorted(allp, reverse=True))),
                    transform=ax.transAxes, fontsize=8.5, color=MUTED, ha='left', va='top')
    legend_maps(fig, 0.075)
    foot(fig, FOOT.replace('\n', ' ').replace(' L = protein length.', '') + '\nLong-range: |i−j| ≥ 24. Chance precision 0.02–0.07. References: crystal structures (PDB) and AlphaFold DB models.')
    fig.savefig(os.path.join(OUT, 'x3%s_matched_%s.png' % (tag, {'a': 'alpha', 'b': 'beta', 'c': 'alphabeta'}[tag])))
    plt.close(fig)

# ---------- 4. summary ----------
ORDER = ['ubiquitin_1UBQ', 'cytb562_256B', 'A0A916UKZ0', 'fn3_tenascin_1TEN', 'A0A8J2YCU1', 'thioredoxin_2TRX', 'A0A5D0P002']
DARK = {'A0A916UKZ0', 'A0A8J2YCU1', 'A0A5D0P002'}
SHORT = {'ubiquitin_1UBQ': 'Ubiquitin', 'cytb562_256B': 'Cyt b562', 'A0A916UKZ0': 'DUF2489',
         'fn3_tenascin_1TEN': 'Tenascin FN3', 'A0A8J2YCU1': 'DUF2958', 'thioredoxin_2TRX': 'Thioredoxin',
         'A0A5D0P002': 'DUF4288'}
xs = [0, 1.6, 2.6, 4.2, 5.2, 6.8, 7.8]
fig = plt.figure(figsize=(8, 4.5), dpi=200)
fig.suptitle('Fable 5.2 contact maps and models, 7 proteins', fontsize=14, color=INK, x=0.03, ha='left', y=0.965)
a1 = fig.add_axes([0.075, 0.29, 0.41, 0.52]); a2 = fig.add_axes([0.575, 0.29, 0.41, 0.52])
for x, t in zip(xs, ORDER):
    rr = [r for r in rows if r['target'] == t and r['job'].startswith('job_')]
    col = NOTPDB if t in DARK else INPDB
    jit = np.linspace(-0.15, 0.15, len(rr)) if len(rr) > 1 else np.zeros(1)
    a1.scatter(x + jit, [float(r['P_sep24_topL']) for r in rr], s=42, c=col, zorder=3, edgecolor='white', lw=1)
    a2.scatter(x + jit, [float(r['TM_best_hand']) for r in rr], s=42, c=col, zorder=3, edgecolor='white', lw=1)
    m = json.load(open(os.path.join(ROOT, 'runs/jobs', rr[0]['job'], 'metrics.json')))
    L = m['L']
    a1.plot([x - 0.35, x + 0.35], [m['native_contacts_sep24'] / ((L - 24) * (L - 23) / 2)] * 2, color=MUTED, lw=2)
    a2.plot([x - 0.35, x + 0.35], [m['TM_ceiling_best_hand']] * 2, color=INK, lw=2)
for ax, ttl in ((a1, 'Long-range precision\n(top-L pairs, |i−j| ≥ 24)'), (a2, 'TM-score of the model built\nfrom Fable\'s contacts')):
    ax.set_xticks(xs); ax.set_xticklabels([SHORT[t] for t in ORDER], fontsize=8.5, rotation=35, ha='right', rotation_mode='anchor')
    ax.set_ylim(0, 1.05); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.tick_params(axis='y', labelsize=9)
    ax.grid(axis='y', color=GRID, lw=0.7); ax.set_axisbelow(True)
    ax.set_title(ttl, fontsize=10.5, color=INK, loc='left')
    for (xa, xb), lab in (((1.6, 2.6), 'all-α'), ((4.2, 5.2), 'all-β'), ((6.8, 7.8), 'α/β')):
        ax.text((xa + xb) / 2, 0.985, lab, ha='center', fontsize=8.5, color=MUTED)
h = [Line2D([], [], marker='o', ls='', color=INPDB, ms=7, label='in PDB, well studied'),
     Line2D([], [], marker='o', ls='', color=NOTPDB, ms=7, label='no relative in PDB'),
     Line2D([], [], color=MUTED, lw=2, label='chance'),
     Line2D([], [], color=INK, lw=2, label='model from the true contacts')]
fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 0.075), ncol=4, fontsize=9, frameon=False)
foot(fig, FOOT.replace(' L = protein length.', '') + ' One dot per run. TM-score: better of the two mirror images.')
fig.savefig(os.path.join(OUT, 'x4_summary.png'))
plt.close(fig)
print('wrote', sorted(os.listdir(OUT)))

# ---------- 5. structures: reference cartoon vs Fable model (CA trace), same superposition ----------
from PIL import Image
SDIR = os.path.join(OUT, 'structures')
sinfo = json.load(open(os.path.join(SDIR, 'structures.json')))


def pair_imgs(t):
    a = Image.open(os.path.join(SDIR, '%s_ref.png' % t)).convert('RGBA')
    b = Image.open(os.path.join(SDIR, '%s_model.png' % t)).convert('RGBA')
    boxes = [im.getchannel('A').getbbox() for im in (a, b)]
    box = (min(x[0] for x in boxes), min(x[1] for x in boxes), max(x[2] for x in boxes), max(x[3] for x in boxes))
    out = []
    for im in (a, b):
        c = im.crop(box)
        bg = Image.new('RGBA', c.size, (255, 255, 255, 255)); bg.alpha_composite(c)
        out.append(np.asarray(bg.convert('RGB')))
    return out


def struct_group(fig, x0, w, y0, h, t, tag, tagcol, name):
    ref, mod = pair_imgs(t)
    reflab = 'AlphaFold DB model' if t.startswith('A0A') else 'crystal structure'
    for k, (img, lab) in enumerate(((ref, reflab), (mod, 'model from Fable\'s contacts'))):
        ax = fig.add_axes([x0 + k * w / 2, y0, w / 2 - 0.01, h])
        ax.imshow(img); ax.axis('off')
        fig.text(x0 + k * w / 2 + (w / 2 - 0.01) / 2, y0 - 0.015, lab, ha='center', va='top', fontsize=9.5, color=MUTED)
    s = sinfo[t]
    fig.text(x0, y0 + h + 0.085, tag, fontsize=10, color=tagcol, fontweight='bold', ha='left')
    fig.text(x0, y0 + h + 0.035, name, fontsize=10.5, color=INK, ha='left')
    note = 'TM-score %.2f' % s['TM']
    if s['n_runs'] > 1:
        note += '  (run with best long-range precision, of %d)' % s['n_runs']
    fig.text(x0 + w / 2, y0 - 0.085, note, fontsize=10, color=INK, ha='center', va='top')


SFOOT = ('Rainbow N→C. Model: Cα positions from the prompt\'s fold() code on Fable\'s final contacts (better mirror image), superposed\n'
         'by TM-score; helices and strands drawn as implied by Fable\'s own contacts, not taken from the reference.')

fig = plt.figure(figsize=(8, 4.5), dpi=200)
fig.suptitle('Ubiquitin (PDB 1UBQ): crystal structure vs model from Fable 5.2 contacts', fontsize=13.5, color=INK, x=0.03, ha='left', y=0.965)
struct_group(fig, 0.2, 0.6, 0.26, 0.6, 'ubiquitin_1UBQ', '', INK, '')
foot(fig, SFOOT)
fig.savefig(os.path.join(OUT, 'x1b_structures_ubiquitin.png'))
plt.close(fig)

for tag, cls, famous, fname, dark, dname in pairs_cls:
    fig = plt.figure(figsize=(8, 4.5), dpi=200)
    fig.suptitle('Same size and class (%s): structures' % cls, fontsize=13.5, color=INK, x=0.03, ha='left', y=0.965)
    struct_group(fig, 0.03, 0.45, 0.26, 0.5, famous, 'in PDB, well studied', INPDB, fname)
    struct_group(fig, 0.52, 0.45, 0.26, 0.5, dark, 'no relative in PDB', NOTPDB, dname)
    foot(fig, SFOOT)
    fig.savefig(os.path.join(OUT, 'x3%s_structures_%s.png' % (tag, {'a': 'alpha', 'b': 'beta', 'c': 'alphabeta'}[tag])))
    plt.close(fig)
print('wrote', sorted(f for f in os.listdir(OUT) if f.endswith('.png')))
