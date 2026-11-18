import numpy as np, json, sys
from lib import *
L = 76
SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

class CS:
    def __init__(self): self.S = {}
    def add(self, i, j, w):
        if i > j: i, j = j, i
        if i < 1 or j > L or j - i < 3: return
        self.S[(i, j)] = max(self.S.get((i, j), 0), w)
    def rm(self, i, j):
        if i > j: i, j = j, i
        self.S.pop((i, j), None)
    def helix(self, a, b, w=0.9):
        for i in range(a, b + 1):
            for k in (3, 4):
                if i + k <= b: self.add(i, i + k, w)
    def anti(self, s, lo, hi, w=1.0, jmin=1, jmax=L, wide=True):
        for i in range(lo, hi + 1):
            j = s - i
            if jmin <= j <= jmax: self.add(i, j, w)
            for dj in (-1, 1):
                if jmin <= j + dj <= jmax: self.add(i, j + dj, w * 0.85)
            if wide:
                for dj in (-2, 2):
                    if jmin <= j + dj <= jmax: self.add(i, j + dj, w * 0.6)
    def par(self, off, lo, hi, w=1.0, jmin=1, jmax=L, wide=True):
        for i in range(lo, hi + 1):
            j = i + off
            if jmin <= j <= jmax: self.add(i, j, w)
            for dj in (-1, 1):
                if jmin <= j + dj <= jmax: self.add(i, j + dj, w * 0.85)
            if wide:
                for dj in (-2, 2):
                    if jmin <= j + dj <= jmax: self.add(i, j + dj, w * 0.6)

def run(cs, tag):
    C = sorted(cs.S)
    print('== %s: %d contacts (%.2f L), %d with sep>=6' % (tag, len(C), len(C) / L, sum(1 for i, j in C if j - i >= 6)))
    X = fold(L, C); r = report(X, C)
    np.save(tag + '.npy', r); np.save(tag + '_X.npy', X)
    cand = [(round(float(r[i - 1, j - 1]), 1), i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1)
            if (i, j) not in cs.S and r[i - 1, j - 1] < 7.0]
    print('close but unpredicted (<7 A):', [(i, j, d) for d, i, j in sorted(cand)[:70]])
    far = sorted([(round(float(r[i - 1, j - 1]), 1), i, j) for (i, j) in C if j - i >= 6], reverse=True)[:15]
    print('longest predicted:', far)
    n = ((r < 8).sum(1) - 1)
    print('neighbours<8A per residue:', ' '.join('%s%d:%d' % (SEQ[k], k + 1, n[k]) for k in range(L)))
    return X, r
