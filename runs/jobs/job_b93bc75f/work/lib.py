import numpy as np

def floyd_warshall(d):
    for m in range(d.shape[0]):
        o = d[m]; d = np.minimum(d, o[None, :] + o[:, None])
    return d

def classic_mds(M):
    L = M.shape[0]; d = np.square(M); c = np.eye(L) - np.ones((L, L)) / L
    s, u = np.linalg.eigh(-0.5 * c @ d @ c)
    return u[:, -3:] * np.sqrt(np.maximum(s[-3:], 0))

def build(L, contacts):
    d = np.full((L, L), np.inf); np.fill_diagonal(d, 0.0)
    k = np.arange(L - 1); d[k, k + 1] = d[k + 1, k] = 3.8
    for i, j in contacts: d[i - 1, j - 1] = d[j - 1, i - 1] = 6.0
    return classic_mds(floyd_warshall(d)), d

def refine(X, contacts, steps=3000, lr=0.01):
    L = len(X); X = X.copy(); B = np.eye(L, k=1, dtype=bool) | np.eye(L, k=-1, dtype=bool)
    C = np.zeros((L, L), bool)
    for i, j in contacts: C[i - 1, j - 1] = C[j - 1, i - 1] = True
    off = ~np.eye(L, dtype=bool)
    for _ in range(steps):
        diff = X[:, None] - X[None]; r = np.linalg.norm(diff, axis=-1) + np.eye(L)
        g = np.where(B, r - 3.8, 0) + np.where(C & (r > 8), r - 8, 0) + np.where(off & ~B & (r < 4), r - 4, 0)
        X -= lr * ((g / r)[..., None] * diff).sum(1)
    return X

def fix_hand(X):
    def dih(p):
        b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]; b1 = b1 / np.linalg.norm(b1)
        v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
        return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))
    h = [dih(X[i:i + 4]) for i in range(len(X) - 3) if np.linalg.norm(X[i] - X[i + 3]) < 6]
    return X * [1, 1, -1] if h and np.median(h) < 0 else X

def fold(L, contacts):
    X, _ = build(L, contacts)
    return fix_hand(refine(fix_hand(X), contacts))

def report(X, contacts):
    L = len(X); r = np.linalg.norm(X[:, None] - X[None], axis=-1)
    bonds = np.diag(r, 1); rg = np.sqrt(((X - X.mean(0)) ** 2).sum(1).mean())
    clash = int(((r < 3.5) & (np.triu(np.ones((L, L)), 3) > 0)).sum())
    bad = sorted([(i, j, round(float(r[i - 1, j - 1]), 1)) for i, j in contacts if r[i - 1, j - 1] > 8], key=lambda t: -t[2])
    print('L %d radius of gyration %.1f A (a compact domain of this length: ~%.1f) bond %.2f +/- %.2f A clashes %d' % (
        L, rg, 2.2 * L ** 0.38, bonds.mean(), bonds.std(), clash))
    print('contacts satisfied (CA-CA <= 8 A): %d of %d' % (len(contacts) - len(bad), len(contacts)))
    print('most violated (i, j, distance): %s' % bad[:30])
    return r

def helix(a, b):
    c = []
    for i in range(a, b + 1):
        for k in (3, 4):
            if i + k <= b: c.append((i, i + k))
    return c

def anti(i0, j0, n, wide=True):
    """antiparallel ladder (i0, j0), (i0+1, j0-1), ... n pairs, plus neighbouring diagonals"""
    c = []
    for k in range(n):
        i, j = i0 + k, j0 - k
        c.append((i, j))
        if wide and k + 1 < n:
            c.append((i + 1, j)); c.append((i, j - 1))
        if wide and k + 2 < n:
            c.append((i, j - 2)); c.append((i + 2, j))
    return c

def par(i0, j0, n):
    c = []
    for k in range(n):
        i, j = i0 + k, j0 + k
        c.append((i, j))
        if k + 1 < n: c.append((i + 1, j)); c.append((i, j + 1))
    return c

def uniq(c):
    return sorted(set((min(a, b), max(a, b)) for a, b in c if abs(a - b) >= 3))

def missed(r, contacts, cut=7.0, sep=6, top=40):
    S = set(contacts); L = len(r); out = []
    for i in range(L):
        for j in range(i + sep, L):
            if r[i, j] < cut and (i + 1, j + 1) not in S: out.append((i + 1, j + 1, round(float(r[i, j]), 1)))
    out.sort(key=lambda t: t[2]); return out[:top]

def segdist(r, segs):
    names = list(segs)
    for a in range(len(names)):
        row = []
        for b in range(len(names)):
            s1, s2 = segs[names[a]], segs[names[b]]
            row.append('%5.1f' % r[s1[0]-1:s1[1], s2[0]-1:s2[1]].min())
        print('%-4s' % names[a], ' '.join(row))
