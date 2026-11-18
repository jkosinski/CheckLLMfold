import json, sys
import numpy as np
from lib import SEQ, fold, report

L = len(SEQ)

def helix_local(a, b):
    c = []
    for i in range(a, b + 1):
        for d in (3, 4):
            if i + d <= b:
                c.append((i, i + d))
    return c

def anti(coreA, flankA, coreB, flankB, S, wide=True):
    c = []
    for i in coreA + flankA:
        for j in coreB + flankB:
            k = abs(i + j - S)
            ic = (i in coreA) + (j in coreB)
            if k <= 1 and ic >= 1:
                c.append((i, j))
            elif wide and k in (3, 4) and ic == 2:
                c.append((i, j))
    return c

def par(coreA, coreB, D):
    return [(i, j) for i in coreA for j in coreB if abs(j - i - D) <= 1]

H = {1: (2, 19), 2: (23, 43), 3: (56, 81), 4: (84, 105)}
core = {1: [3, 7, 10, 14, 17], 2: [23, 26, 30, 33, 37, 40], 3: [58, 61, 65, 68, 72, 75, 79], 4: [84, 87, 91, 94, 98, 101, 105]}
flank = {1: [6, 13, 16, 20], 2: [29, 36, 38, 34, 27], 3: [62, 69, 76, 64], 4: [90, 97, 88, 95, 102]}

def make(round_no):
    c = []
    for h in H.values():
        c += helix_local(*h)
    # helix 1 - helix 2 (antiparallel, i + j ~ 43/44)
    c += anti(core[1], flank[1], core[2], flank[2], 43.5)
    # helix 3 - helix 4 (antiparallel, i + j ~ 166)
    c += anti(core[3], flank[3], core[4], flank[4], 166)
    # turns
    c += [(19, 22), (20, 23), (19, 23), (20, 24), (21, 24), (80, 83), (81, 84), (80, 84), (81, 85), (79, 84), (18, 24), (78, 85)]
    if round_no >= 1:
        # helix 2 - helix 3 (antiparallel, i + j ~ 98)
        c += anti(core[2], flank[2], core[3], flank[3], 98, wide=False)
        # helix 1 - helix 4 (antiparallel, i + j ~ 108)
        c += anti(core[1], flank[1], core[4], flank[4], 108, wide=False)
    return sorted(set((min(a, b), max(a, b)) for a, b in c if abs(a - b) >= 3))

if __name__ == '__main__':
    c = make(1)
    print(len(c), 'contacts; target', 2 * L, '-', 3 * L)
    X = fold(L, c)
    r = report(X, c)
    np.save('X1.npy', X)
    json.dump(c, open('c1.json', 'w'))
    S = set(c)
    cand = [(round(float(r[i - 1, j - 1]), 1), i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1) if (i, j) not in S and r[i - 1, j - 1] < 7.5]
    print('close but unpredicted (|i-j|>=6):', sorted(cand)[:80])
    # helix-helix centroid distances
    for a in H:
        for b in H:
            if a < b:
                ca = X[H[a][0] - 1:H[a][1]].mean(0); cb = X[H[b][0] - 1:H[b][1]].mean(0)
                print('helix', a, b, 'centroid dist %.1f' % np.linalg.norm(ca - cb))
