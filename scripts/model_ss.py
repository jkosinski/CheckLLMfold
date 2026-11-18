"""Secondary structure implied by Fable's final contact list, for drawing its CA model as a ribbon.
Strand: residue pairs (i, j), j - i >= 6, that sit in a ladder of predicted contacts: (i-1, j+1) and (i+1, j-1)
(antiparallel) or (i-1, j-1) and (i+1, j+1) (parallel); runs of length >= 2 kept.
Helix: three consecutive helical turns, (k+n, k+n+3) and (k+n, k+n+4) for n = 0..2, outside strands; runs >= 6 kept.
Writes results/x_thread/structures/<target>_model_ss.json and prints the SS strings next to the reference DSSP-like SS."""
import os, sys, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from evaluate import last_json
from structutil import parse
from sscomp import dssp_lite

SDIR = os.path.join(ROOT, 'results/x_thread/structures')
info = json.load(open(os.path.join(SDIR, 'structures.json')))
REFS = {'ubiquitin_1UBQ': ('data/1UBQ/1ubq.pdb', 'A'), 'cytb562_256B': ('data/256B/256b.pdb', 'A'),
        'fn3_tenascin_1TEN': ('data/1TEN/1ten.pdb', 'A'), 'thioredoxin_2TRX': ('data/2TRX/2trx.pdb', 'A')}


def ss_from_contacts(L, contacts):
    C = set(tuple(sorted(c)) for c in contacts)
    has = lambda i, j: (min(i, j), max(i, j)) in C
    ss = ['-'] * (L + 2)
    # strands first: ladders of predicted contacts with |i - j| >= 6
    for (i, j) in C:
        if j - i < 6:
            continue
        if (has(i - 1, j + 1) and has(i + 1, j - 1)) or (has(i - 1, j - 1) and has(i + 1, j + 1)):
            ss[i] = ss[j] = 'E'
    # helices: three consecutive helical turns (k+n, k+n+3) and (k+n, k+n+4), n = 0..2, outside strands
    hel = [False] * (L + 2)
    for k in range(1, L - 5):
        if all(has(k + n, k + n + 3) and has(k + n, k + n + 4) for n in range(3)):
            for m in range(k, k + 7):
                hel[m] = True
    for m in range(1, L + 1):
        if hel[m] and ss[m] != 'E':
            ss[m] = 'H'
    # break strand runs at hairpin turns: an antiparallel ladder pair (i, j) with both ends in one strand run
    runs, k = [], 1
    while k <= L:
        if ss[k] == 'E':
            m = k
            while m + 1 <= L and ss[m + 1] == 'E':
                m += 1
            runs.append((k, m)); k = m + 1
        else:
            k += 1
    for a, b in runs:
        pairs = [(i, j) for (i, j) in C if a <= i and j <= b and j - i >= 5 and
                 any(has(i + d, j - d) for d in (-2, -1, 1, 2))]
        if pairs:
            i, j = min(pairs, key=lambda p: p[1] - p[0])
            c = (i + j) // 2
            for m in (c, c + 1) if (i + j) % 2 else (c - 1, c, c + 1):
                ss[m] = '-'
    s = ''.join(ss[1:L + 1])
    # drop helix runs shorter than 6 and strand runs shorter than 2
    out = list(s)
    for ch, minlen in (('H', 6), ('E', 2)):
        k = 0
        while k < len(out):
            if out[k] == ch:
                m = k
                while m + 1 < len(out) and out[m + 1] == ch:
                    m += 1
                if m - k + 1 < minlen:
                    for q in range(k, m + 1):
                        out[q] = '-'
                k = m + 1
            else:
                k += 1
    return ''.join(out)


def segments(s, ch):
    out, k = [], 0
    while k < len(s):
        if s[k] == ch:
            m = k
            while m + 1 < len(s) and s[m + 1] == ch:
                m += 1
            out.append((k + 1, m + 1)); k = m + 1
        else:
            k += 1
    return out


for t, d in info.items():
    ans = last_json(open(os.path.join(ROOT, 'runs/jobs', d['job'], 'final.txt')).read())
    s = ss_from_contacts(d['L'], ans['final_contacts'])
    f, ch = REFS.get(t, ('data/%s/%s.pdb' % (t, t), None))
    ref = dssp_lite([r for r in parse(os.path.join(ROOT, f), ch) if 'CA' in r['atoms']])
    json.dump({'ss': s, 'helix': segments(s, 'H'), 'strand': segments(s, 'E')}, open(os.path.join(SDIR, '%s_model_ss.json' % t), 'w'))
    print('%-18s ref   %s\n%-18s fable %s' % (t, ref, '', s))
