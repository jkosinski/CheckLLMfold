"""Figures (1600x900) for an X thread on next-amino-acid prediction by Claude Fable 5.2. Same style as
../scripts/figures_x.py. Run with the data_analysis conda env (matplotlib). Writes results/x_thread/y1..y4 PNGs."""
import os, csv, json, gzip, math, random, collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'results/x_thread'); os.makedirs(OUT, exist_ok=True)
INK, MUTED, GRID, BG = '#0b0b0b', '#5f5e5a', '#e6e5e0', '#ffffff'
OLD, NEW = '#6250d6', '#1baf7a'          # public before the training cutoff / first public after it
NEWRES, OLDRES, OTHER = '#2a78d6', '#eb6834', '#cfcec8'
plt.rcParams.update({'font.size': 13, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': BG, 'savefig.facecolor': BG, 'font.family': 'DejaVu Sans'})
W, H = 16, 9


def load(path):
    op = gzip.open if path.endswith('.gz') else open
    return [json.loads(l) for l in op(os.path.join(ROOT, path), 'rt')]


def wilson(k, n):
    if not n: return (0, 0)
    p = k / n; z = 1.96; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


def fam_gap(rows):
    fam = collections.defaultdict(lambda: {'real': [], 'shuffled': []})
    for o in rows: fam[o['pfam']][o['condition']].append(bool(o['first_correct']))
    g = [sum(v['real']) / len(v['real']) - sum(v['shuffled']) / len(v['shuffled']) for v in fam.values() if v['real'] and v['shuffled']]
    m = sum(g) / len(g); rnd = random.Random(0)
    bs = sorted(sum(rnd.choice(g) for _ in g) / len(g) for _ in range(4000))
    return 100 * m, 100 * bs[100], 100 * bs[3899]


def foot(fig, text):
    fig.text(0.012, 0.012, text, fontsize=10, color=MUTED, ha='left', va='bottom')


def title(fig, t, sub=None):
    fig.text(0.03, 0.93, t, fontsize=24, color=INK, ha='left', va='bottom')
    if sub: fig.text(0.03, 0.885, sub, fontsize=14, color=MUTED, ha='left', va='bottom')


# ---------- y1: well-known proteins ----------
def y1():
    rows = [o for o in load('results/famous_fable52.jsonl') if 'error' not in o]
    names = {'P0CG48': 'Ubiquitin', 'P42212': 'GFP', 'P00698': 'Lysozyme (hen)', 'P68871': 'Hemoglobin β',
             'P02185': 'Myoglobin (whale)', 'P01116': 'KRAS', 'P04637': 'p53', 'P01308': 'Insulin',
             'P0AA25': 'Thioredoxin (E. coli)', 'P00441': 'SOD1', 'P00974': 'BPTI', 'P0DP23': 'Calmodulin'}
    order = list(names)
    fig = plt.figure(figsize=(W, H), dpi=100)
    title(fig, 'Well-known proteins: Fable 5.2 recites them',
          'Next residue after 5 random cut points per protein, real sequence vs the same residues shuffled')
    ax = fig.add_axes([0.2, 0.17, 0.5, 0.66])
    for y, acc in enumerate(order):
        for cond, col, mk in (('real', OLD, 'o'), ('shuffled', MUTED, 'o')):
            rs = [o for o in rows if o['accession'] == acc and o['condition'] == cond]
            v = 100 * sum(o['first_correct'] for o in rs) / len(rs)
            ax.scatter(v, y, s=140, color=col if cond == 'real' else BG, edgecolor=col, linewidth=2, zorder=3)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([names[a] for a in order], color=INK)
    ax.set_ylim(len(order) - 0.5, -0.5); ax.set_xlim(-5, 105)
    ax.set_xlabel('next residue right (% of 5 cut points)')
    ax.grid(axis='x', color=GRID); ax.set_axisbelow(True)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=11, color=OLD, label='real sequence'),
                       Line2D([], [], marker='o', ls='', ms=11, mfc=BG, mec=MUTED, mew=2, label='shuffled')],
              loc='lower center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
    R = [o for o in rows if o['condition'] == 'real']; S = [o for o in rows if o['condition'] == 'shuffled']
    x0 = 0.74
    for k, (big, small) in enumerate([('%d / %d' % (sum(o['first_correct'] for o in R), len(R)), 'next residue right'),
                                      ('%d / %d' % (sum(o['identity'] == 1 for o in R), len(R)), 'all 10 next residues\nexactly right'),
                                      ('%.0f%%' % (100 * sum(o['first_correct'] for o in S) / len(S)), 'next residue right\non shuffled sequences')]):
        y = 0.70 - k * 0.2
        fig.text(x0, y, big, fontsize=34, color=INK, ha='left', va='bottom')
        fig.text(x0, y - 0.01, small, fontsize=13, color=MUTED, ha='left', va='top')
    foot(fig, 'Claude Fable 5.2 (effort low) through headless Claude Code, no tools. Prompt as in Tu & Chen 2026: '
              '"Reply with exactly the next 10 amino acids".')
    fig.savefig(os.path.join(OUT, 'y1_famous.png')); plt.close(fig)


# ---------- y2: gap and exact runs by number of close relatives (Tu & Chen's published replies) ----------
def y2():
    ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open(os.path.join(ROOT, 'data/main_uniref50.tsv')), delimiter='\t')}
    def stratum(a):
        c = ur50.get(a)
        if c is None: return 4
        return 3 if c == 1 else (2 if c <= 5 else (1 if c <= 50 else 0))
    labels = ['50+', '5–49', '1–4', 'none\npublic before\ncutoff', 'none\nnew after\ncutoff']
    cols = [OLD, OLD, OLD, OLD, NEW]
    models = [('Claude Fable 5.2', 'data/tu_outputs/confirm/claude-fable-5-2.jsonl.gz', 'data/tu_outputs/d_set/claude-fable-5-2.jsonl.gz'),
              ('ProGen3-3B', 'data/tu_outputs/plm_confirm/progen3-3b__prefix.jsonl.gz', 'data/tu_outputs/plm_d/progen3-3b__prefix.jsonl.gz')]
    fig = plt.figure(figsize=(W, H), dpi=100)
    title(fig, 'The fewer known relatives, the worse Fable 5.2 does',
          '668 + 56 little-studied bacterial proteins, Tu & Chen\'s published replies, split by close relatives in UniProt')
    ax = fig.add_axes([0.07, 0.24, 0.4, 0.55]); bx = fig.add_axes([0.57, 0.24, 0.4, 0.55])
    for mi, (name, a, b) in enumerate(models):
        rows = [o for o in load(a) + load(b) if o.get('valid', True) is not False]
        for s in range(5):
            m, lo, hi = fam_gap([o for o in rows if stratum(o['accession']) == s])
            x = s + (-0.12 if mi == 0 else 0.12)
            ax.plot([x, x], [lo, hi], color=cols[s] if mi == 0 else MUTED, lw=2, zorder=2)
            ax.scatter(x, m, s=150, color=cols[s] if mi == 0 else BG, edgecolor=cols[s] if mi == 0 else MUTED, lw=2, zorder=3)
            if mi == 0: ax.text(x - 0.1, m, '%.0f' % m, ha='right', va='center', fontsize=12, color=INK)
        if mi == 0:
            for s in range(5):
                rs = [o for o in rows if stratum(o['accession']) == s and o['condition'] == 'real']
                k = sum(o['exact_run'] >= 5 for o in rs); n = len(rs); lo, hi = wilson(k, n)
                bx.bar(s, 100 * k / n, width=0.62, color=cols[s], zorder=2)
                bx.plot([s, s], [lo, hi], color=INK, lw=1.2, zorder=3)
                bx.text(s, hi + 0.6, '%.1f%%' % (100 * k / n), ha='center', va='bottom', fontsize=12, color=INK)
    for z in (ax, bx):
        z.set_xticks(range(5)); z.set_xticklabels(labels, fontsize=12, color=INK)
        z.set_xlabel('relatives in UniProt at ≥50% identity', labelpad=8)
        z.grid(axis='y', color=GRID); z.set_axisbelow(True); z.set_xlim(-0.6, 4.6)
    ax.set_ylim(0, 60); bx.set_ylim(0, 22)
    ax.set_title('Next residue: real minus shuffled accuracy (points)', loc='left', fontsize=14, color=INK, pad=12)
    bx.set_title('Cut points where the first 5 residues are exactly right (%)', loc='left', fontsize=14, color=INK, pad=12)
    ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=11, color=OLD, label='Fable 5.2, public before cutoff'),
                       Line2D([], [], marker='o', ls='', ms=11, color=NEW, label='Fable 5.2, new after cutoff'),
                       Line2D([], [], marker='o', ls='', ms=11, mfc=BG, mec=MUTED, mew=2, label='ProGen3-3B (protein language model)')],
              loc='upper right', fontsize=11, frameon=False)
    foot(fig, 'Data and replies: github.com/XinmingTu/LLM-bioseq-probe. "New" = 56 RefSeq proteins first public in August 2026, after '
              'Fable 5.2\'s June 2026 training cutoff, with no earlier\nrelative at ≥50% identity. Gap: mean over Pfam families, 95% '
              'bootstrap CI. Relatives = other members of the protein\'s UniRef50 cluster.')
    fig.savefig(os.path.join(OUT, 'y2_relatives.png')); plt.close(fig)


# ---------- y3: matched old vs new, both setups ----------
def y3():
    ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open(os.path.join(ROOT, 'data/main_uniref50.tsv')), delimiter='\t')}
    ours = [o for o in load('results/matched_fable52.jsonl') if 'error' not in o and o.get('valid')]
    tu = [o for o in load('data/tu_outputs/confirm/claude-fable-5-2.jsonl.gz') if o.get('valid') and ur50[o['accession']] == 1]
    tu_d = [o for o in load('data/tu_outputs/d_set/claude-fable-5-2.jsonl.gz') if o.get('valid')]
    groups = [('Tu & Chen\'s runs\n(API)', tu, tu_d), ('This test\n(headless Claude Code)', [o for o in ours if o['set'] == 'pre_single'],
                                                       [o for o in ours if o['set'] == 'post_D'])]
    fig = plt.figure(figsize=(W, H), dpi=100)
    title(fig, 'Same distance to known relatives, different outcome',
          'Proteins with no relative at ≥50% identity: public before Fable 5.2\'s cutoff vs first public after it')
    ax = fig.add_axes([0.25, 0.2, 0.62, 0.6])
    ys = []
    for gi, (lab, old, new) in enumerate(groups):
        for k, (rows, col, nm) in enumerate(((old, OLD, 'public before cutoff'), (new, NEW, 'new after cutoff'))):
            y = gi * 3 + k; ys.append((y, nm))
            m, lo, hi = fam_gap(rows)
            ax.plot([lo, hi], [y, y], color=col, lw=2.5); ax.scatter(m, y, s=170, color=col, zorder=3)
            ax.text(hi + 1, y, '%.1f' % m, va='center', fontsize=13, color=INK)
        ax.annotate(lab, xy=(0.03, gi * 3 + 0.5), xycoords=('figure fraction', 'data'), fontsize=14, color=INK, va='center')
    ax.set_yticks([y for y, _ in ys]); ax.set_yticklabels([n for _, n in ys], color=INK)
    ax.set_ylim(4.6, -0.6); ax.set_xlim(0, 50)
    ax.set_xlabel('next residue: real minus shuffled accuracy (points), 95% CI')
    ax.grid(axis='x', color=GRID); ax.set_axisbelow(True)
    foot(fig, 'Before cutoff: the 116 proteins of Tu & Chen\'s main set alone in their UniRef50 cluster. After cutoff: their 56 RefSeq '
              'proteins first public in August 2026 (closest earlier relative 34-49%).\nThis test asked 2 cut points per old protein '
              'and 5 per new one, shuffled controls on 1-2. Headless Claude Code scores lower than the API overall.')
    fig.savefig(os.path.join(OUT, 'y3_matched.png')); plt.close(fig)


# ---------- y4: divergence test ----------
def y4():
    d = os.path.join(ROOT, 'data/div_uniprot')
    qs = list(csv.DictReader(open(os.path.join(d, 'q_div.tsv')), delimiter='\t'))
    rq = [r for r in csv.DictReader(open(os.path.join(d, 'q_rel.tsv')), delimiter='\t') if r['condition'] == 'real']
    new = {o['qid']: o for o in load('results/divU_fable52.jsonl')}
    rel = {(o['new_accession'], o['site'], int(o['cut'])): o for o in load('results/relselfU_fable52.jsonl') if o['condition'] == 'real'}
    P = [(a, new[a['qid']], rel[(b['new_accession'], b['site'], int(b['cut']))]) for a, b in zip(qs, rq)]
    fig = plt.figure(figsize=(W, H), dpi=100)
    title(fig, 'On a new protein, Fable 5.2 gives the residue of the old one',
          '48 proteins first public after Fable 5.2\'s cutoff, each 70-95% identical to an older UniProt protein')
    # example
    mono = {'family': 'DejaVu Sans Mono', 'fontsize': 18}
    ex = [('older UniProt protein (A0A7X1GB55, 2021)', '…HEVSKRQADLDKAMKKGDPE', 'K', 'INKRKDKLA', OLDRES),
          ('new protein (WP_526858458, first public 12 Nov 2026)', '…HEVSERQADLDKAMKKGDPE', 'R', 'INKRKDKLA', NEWRES),
          ('Fable 5.2, given the new protein up to the cut', '', 'K', 'INKRKDKLA', OLDRES)]
    for k, (lab, pre, nxt, rest, col) in enumerate(ex):
        y = 0.775 - k * 0.07
        fig.text(0.03, y, lab, fontsize=13, color=MUTED, va='center')
        fig.text(0.665, y, pre, color=INK, va='center', ha='right', **mono)
        fig.text(0.672, y, nxt, color=col, va='center', ha='left', fontweight='bold', **mono)
        fig.text(0.69, y, rest, color=INK, va='center', ha='left', **mono)
    fig.add_artist(Line2D([0.6685, 0.6685], [0.62, 0.815], color=MUTED, lw=1, ls=':'))
    fig.text(0.6685, 0.825, 'cut', fontsize=11, color=MUTED, ha='center', va='bottom')
    # bars
    ax = fig.add_axes([0.3, 0.2, 0.55, 0.26])
    groups = [('Fable knows the old residue', lambda r: r['first_correct']), ('Fable does not', lambda r: not r['first_correct'])]
    for y, (lab, sel) in enumerate(groups):
        S = [(a, n, r) for a, n, r in P if a['site'] == 'div' and sel(r)]
        m = len(S)
        fn = 100 * sum(n['pred'][0] == a['true_res'] for a, n, r in S) / m
        fo = 100 * sum(n['pred'][0] == a['rel_res'] for a, n, r in S) / m
        x = 0
        for v, col, txtcol in ((fn, NEWRES, BG), (fo, OLDRES, INK), (100 - fn - fo, OTHER, INK)):
            ax.barh(y, v - 0.4, left=x, height=0.62, color=col, zorder=2)
            if v >= 6: ax.text(x + (v - 0.4) / 2, y, '%.0f%%' % v, ha='center', va='center', fontsize=13, color=txtcol)
            x += v
        ax.text(-1.5, y, '%s (n=%d)' % (lab, m), ha='right', va='center', fontsize=13, color=INK)
    ax.set_xlim(0, 100); ax.set_ylim(1.6, -0.6); ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('what Fable gives at a position where the new protein changed (%)')
    ax.text(100, -0.62, '"knows" = given the old protein, Fable predicts its residue here', ha='right', fontsize=11, color=MUTED, transform=ax.transData)
    ax.legend(handles=[Patch(color=NEWRES, label='new protein\'s residue'), Patch(color=OLDRES, label='old protein\'s residue'),
                       Patch(color=OTHER, label='something else')], loc='lower left', bbox_to_anchor=(-0.02, 1.12), ncol=3,
              frameon=False, fontsize=12)
    cons = [(a, n, r) for a, n, r in P if a['site'] == 'cons']
    k1 = [(a, n, r) for a, n, r in cons if r['first_correct']]
    fig.text(0.03, 0.085, 'Where nothing changed, Fable is right on the new protein at %d%% of positions where it knows the old residue '
             '(n=%d).' % (round(100 * sum(n['pred'][0] == a['true_res'] for a, n, r in k1) / len(k1)), len(k1)), fontsize=12, color=INK)
    foot(fig, 'Claude Fable 5.2 (effort low), headless Claude Code. New proteins: RefSeq, created Aug-Nov 2026, exact sequence not in UniParc before July 2026.\n'
              'Old relatives: UniProtKB entries created before June 2026, found with MMseqs2. 3 changed and 1 unchanged position per protein.')
    fig.savefig(os.path.join(OUT, 'y4_divergence.png')); plt.close(fig)

# ---------- y3b: seen vs new, Fable against protein language models ----------
def y3b():
    ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open(os.path.join(ROOT, 'data/main_uniref50.tsv')), delimiter='\t')}
    T = 'data/tu_outputs/'
    def tu(a, b):
        old = [o for o in load(T + a + '.jsonl.gz') if o.get('valid', True) is not False and ur50[o['accession']] == 1]
        new = [o for o in load(T + b + '.jsonl.gz') if o.get('valid', True) is not False]
        return old, new
    ours = [o for o in load('results/matched_fable52.jsonl') if 'error' not in o and o.get('valid')]
    rows = [('Claude Fable 5.2', *tu('confirm/claude-fable-5-2', 'd_set/claude-fable-5-2')),
            ('Claude Opus 5.5', *tu('confirm/claude-opus-5-5', 'd_set/claude-opus-5-5')),
            ('ProGen3-3B', *tu('plm_confirm/progen3-3b__prefix', 'plm_d/progen3-3b__prefix')),
            ('ProtGPT3-10B', *tu('plm_confirm/protgpt3-10B__prefix', 'plm_d/protgpt3-10B__prefix')),
            ('Claude Fable 5.2, this test', [o for o in ours if o['set'] == 'pre_single'], [o for o in ours if o['set'] == 'post_D'])]
    fig = plt.figure(figsize=(W, H), dpi=100)
    title(fig, 'On unseen proteins Fable 5.2 keeps 61% of its gap, protein language models 35-51%',
          'Proteins with no relative at ≥50% identity, public before Fable 5.2\'s cutoff vs first public after it')
    ax = fig.add_axes([0.25, 0.2, 0.55, 0.6])
    for y, (name, old, new) in enumerate(rows):
        yy = y + (0.5 if y == 4 else 0)
        mo, lo_o, hi_o = fam_gap(old); mn, lo_n, hi_n = fam_gap(new)
        ax.plot([mn, mo], [yy, yy], color=GRID, lw=6, solid_capstyle='round', zorder=1)
        for m, lo, hi, col, dy in ((mo, lo_o, hi_o, OLD, -0.13), (mn, lo_n, hi_n, NEW, 0.13)):
            ax.plot([lo, hi], [yy + dy, yy + dy], color=col, lw=1.5, alpha=0.7, zorder=2)
            ax.scatter(m, yy + dy, s=150, color=col, zorder=3)
        ax.text(52, yy, '%.0f%%' % (100 * mn / mo), va='center', ha='left', fontsize=14, color=INK)
        ax.text(-1.5, yy, name, va='center', ha='right', fontsize=14, color=INK)
    ax.text(52, -0.85, 'new as % of old', ha='left', va='center', fontsize=12, color=MUTED)
    ax.set_yticks([]); ax.spines['left'].set_visible(False)
    ax.set_ylim(5.0, -1.1); ax.set_xlim(0, 50)
    ax.set_xlabel('next residue: real minus shuffled accuracy (points), 95% CI')
    ax.grid(axis='x', color=GRID); ax.set_axisbelow(True)
    ax.axhline(3.75, color=MUTED, lw=0.6, ls=':')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=11, color=OLD, label='public before cutoff (116 proteins)'),
                       Line2D([], [], marker='o', ls='', ms=11, color=NEW, label='first public after cutoff (56 proteins)')],
              loc='lower left', bbox_to_anchor=(-0.02, 1.02), ncol=2, frameon=False, fontsize=12)
    foot(fig, 'Top four rows: Tu & Chen\'s published replies (github.com/XinmingTu/LLM-bioseq-probe). 85 of the 116 older proteins were in UniProt by April 2021,\n'
              'so they are likely in the protein language models\' training data too. Bottom row: this test, headless Claude Code, which scores lower than the API.')
    fig.savefig(os.path.join(OUT, 'y3b_vs_plm.png')); plt.close(fig)


if __name__ == '__main__':
    y1(); y2(); y3(); y3b(); y4()
    print('written to', OUT)
