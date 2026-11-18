"""Crude DSSP (Kabsch-Sander H-bond energy) for backbone N, CA, C, O. Returns string of H/E/-."""
import numpy as np


def dssp_lite(res):
    N = np.array([r['atoms']['N'] for r in res]); CA = np.array([r['atoms']['CA'] for r in res])
    C = np.array([r['atoms']['C'] for r in res]); O = np.array([r['atoms']['O'] for r in res])
    L = len(res)
    H = N.copy()
    for i in range(1, L):
        v = C[i - 1] - O[i - 1]; H[i] = N[i] + v / np.linalg.norm(v)
    # E(i->j): C=O of i accepts from N-H of j
    def dist(a, b): return np.linalg.norm(a[:, None] - b[None], axis=-1)
    rON, rCH, rOH, rCN = dist(O, N), dist(C, H), dist(O, H), dist(C, N)
    E = 0.084 * 332 * (1 / rON + 1 / rCH - 1 / rOH - 1 / rCN)
    hb = E < -0.5
    np.fill_diagonal(hb, False)
    for i in range(L):
        hb[i, max(0, i - 1):i + 2] = False
    hb[:, 0] = False  # residue 1 has no H
    ss = ['-'] * L
    # helices: two consecutive 4-turns
    turn4 = [i + 4 < L and hb[i, i + 4] for i in range(L)]
    for i in range(1, L - 4):
        if turn4[i - 1] and turn4[i]:
            for k in range(i, i + 4):
                ss[k] = 'H'
    # bridges
    def hbond(a, b): return 0 <= a < L and 0 <= b < L and hb[a, b]
    bridge = [False] * L
    for i in range(1, L - 1):
        for j in range(1, L - 1):
            if abs(i - j) < 3:
                continue
            par = (hbond(i - 1, j) and hbond(j, i + 1)) or (hbond(j - 1, i) and hbond(i, j + 1))
            anti = (hbond(i, j) and hbond(j, i)) or (hbond(i - 1, j + 1) and hbond(j - 1, i + 1))
            if par or anti:
                bridge[i] = True
    for i in range(L):
        if bridge[i] and ss[i] != 'H':
            ss[i] = 'E'
    # fill single gaps between E (ladders with bulges)
    for i in range(1, L - 1):
        if ss[i] == '-' and ss[i - 1] == 'E' and ss[i + 1] == 'E':
            ss[i] = 'E'
    return ''.join(ss)
