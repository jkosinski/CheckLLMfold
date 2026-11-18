import numpy as np, json, sys
from lib import *

L = 107
seq = "MGRAGWYGVRCVFRWVHEGRQVYEETVTVWRAESFEEALQKAEVGAFEYAAECDGQYLEFAQAYFIGEDEVIGEGAEVFSLMRESELGGRDYVTRYFDTGDERQGAW"
assert len(seq) == L


def helix(a, b):
    return [(i, i + k) for i in range(a, b + 1) for k in (3, 4) if i + k <= b]


def anti(s, lo, hi, w=1):
    """antiparallel pairing i + j = s for i in lo..hi, plus register neighbours j +/- w"""
    out = []
    for i in range(lo, hi + 1):
        j = s - i
        for dj in range(-w, w + 1):
            if j + dj - i >= 3:
                out.append((i, j + dj))
    return out


def par(off, lo, hi, w=1):
    out = []
    for i in range(lo, hi + 1):
        for dj in range(-w, w + 1):
            out.append((i, i + off + dj))
    return out


def clean(c):
    s = set()
    for i, j in c:
        i, j = min(i, j), max(i, j)
        if j - i >= 3 and 1 <= i and j <= L:
            s.add((i, j))
    return sorted(s)


def run(c, name):
    c = clean(c)
    print(name, 'n contacts', len(c), 'long(>=6)', sum(1 for i, j in c if j - i >= 6))
    X = fold(L, c)
    r = report(X, c)
    np.save(name + '_r.npy', r)
    json.dump(c, open(name + '_c.json', 'w'))
    cs = set(c)
    miss = sorted([(round(float(r[i - 1, j - 1]), 1), i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1)
                   if (i, j) not in cs and r[i - 1, j - 1] < 7.0])
    print('unpredicted close pairs (<7A, sep>=6):', len(miss))
    print([(i, j, d) for d, i, j in miss[:80]])
    # per-residue neighbour count at 8 A (sep >= 3)
    nb = [(int(((r[i] < 8).sum() - 1 - sum(1 for k in (-2, -1, 1, 2) if 0 <= i + k < L)))) for i in range(L)]
    print('neighbour counts:', ' '.join('%d%s:%d' % (i + 1, seq[i], nb[i]) for i in range(L)))
    return X, r, c
