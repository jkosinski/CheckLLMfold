import json, sys
import numpy as np
from lib import SEQ, fold, report
from r1 import helix_local, anti, par, H, core, flank

L = len(SEQ)

def norm(c):
    return sorted(set((min(a, b), max(a, b)) for a, b in c if abs(a - b) >= 3))

def make(n):
    c = []
    for h in H.values():
        c += helix_local(*h)
    c += anti(core[1], flank[1], core[2], flank[2], 43.5)
    c += anti(core[3], flank[3], core[4], flank[4], 166)
    c += [(19, 22), (20, 23), (19, 23), (20, 24), (21, 24), (80, 83), (81, 84), (80, 84), (81, 85), (79, 84), (18, 24), (78, 85)]
    c += anti(core[2], flank[2], core[3], flank[3], 98, wide=(n >= 2))
    c += anti(core[1], flank[1], core[4], flank[4], 108, wide=(n >= 2))
    if n >= 2:
        # core-core pairs between helix 1 / helix 3 (same direction, j - i ~ 55) and helix 2 / helix 4 (j - i ~ 64.5)
        c += par(core[1], core[3], 55)
        c += par(core[2], core[4], 64.5)
        # loop 43-56: local turns and L48 packing against end of helix 2 / start of helix 3
        c += [(42, 45), (43, 46), (44, 47), (45, 48), (46, 49), (47, 50), (48, 51), (49, 52), (50, 53), (51, 54), (52, 55),
              (53, 56), (54, 57), (55, 58), (54, 58), (55, 59), (53, 57),
              (48, 38), (48, 41), (48, 58), (48, 61), (44, 58), (41, 58), (44, 55), (41, 55), (45, 54)]
    if n >= 3:
        # model put these close in round 2 and they are chemically sensible (apolar / same face of the helices)
        c += EXTRA3
    if n >= 4:
        c += EXTRA4
        c[:] = [p for p in norm(c) if p not in set(norm(DROP4))]
    return norm(c)

EXTRA3 = []
EXTRA4 = []
DROP4 = []

def run(n, c=None):
    c = c or make(n)
    print('round', n, ':', len(c), 'contacts; allowed', 2 * L, '-', 3 * L)
    X = fold(L, c)
    r = report(X, c)
    np.save('X%d.npy' % n, X); np.save('r%d.npy' % n, r)
    json.dump(c, open('c%d.json' % n, 'w'))
    S = set(c)
    cand = [(round(float(r[i - 1, j - 1]), 1), i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1) if (i, j) not in S and r[i - 1, j - 1] < 7.0]
    print('close but unpredicted (|i-j|>=6):', len(cand), sorted(cand)[:120])
    far = [(round(float(r[i - 1, j - 1]), 1), i, j) for i, j in c if abs(i - j) >= 6 and r[i - 1, j - 1] > 7.5]
    print('predicted, strained (>7.5):', sorted(far, reverse=True)[:60])
    for a in H:
        for b in H:
            if a < b:
                ca = X[H[a][0] - 1:H[a][1]].mean(0); cb = X[H[b][0] - 1:H[b][1]].mean(0)
                print('helix', a, b, 'centroid dist %.1f' % np.linalg.norm(ca - cb))
    return c, X, r

if __name__ == '__main__':
    run(int(sys.argv[1]))
