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

SEQ = 'SDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA'
L = len(SEQ)

def helix(a, b):
    s = set()
    for i in range(a, b + 1):
        for k in (3, 4):
            if i + k <= b: s.add((i, i + k))
    return s

def par(a, b, off, w=1):
    s = set()
    for i in range(a, b + 1):
        j = i + off
        s.add((i, j))
        if w: s.add((i, j + 1)); s.add((i, j - 1))
    return s

def anti(a, b, sm, w=1):
    s = set()
    for i in range(a, b + 1):
        j = sm - i
        if j - i < 3: continue
        s.add((i, j))
        if w:
            if j - 1 - i >= 3: s.add((i, j - 1))
            s.add((i, j + 1))
    return s

def norm(s):
    out = set()
    for i, j in s:
        if i > j: i, j = j, i
        if j - i >= 3 and 1 <= i and j <= L: out.add((i, j))
    return sorted(out)

def unpred(r, contacts, cut=7.0, minsep=6, n=45):
    C = set(contacts)
    c = [(round(float(r[i-1, j-1]), 1), i, j) for i in range(1, L+1) for j in range(i+minsep, L+1)
         if r[i-1, j-1] < cut and (i, j) not in C]
    c.sort()
    print('close but unpredicted (%d):' % len(c), [(i, j, d) for d, i, j in c[:n]])

SEGS = {'b1': (3, 8), 'a1': (11, 19), 'b2': (23, 28), 'a2': (34, 49), 'b3': (53, 59),
        'a3': (63, 70), 'b4': (76, 82), 'b5': (85, 91), 'a4': (95, 107)}

def segdist(r):
    names = list(SEGS)
    print('min distance between segments:')
    print('     ' + ' '.join('%5s' % n for n in names))
    for a in names:
        row = []
        for b in names:
            (a1, a2), (b1, b2) = SEGS[a], SEGS[b]
            row.append('%5.1f' % r[a1-1:a2, b1-1:b2].min() if a != b else '   - ')
        print('%4s ' % a + ' '.join(row))

def run(s, tag):
    c = norm(s)
    print(tag, 'L', L, '|', len(c), 'contacts (%.2f x L); |i-j|>=6: %d' % (len(c) / L, sum(1 for i, j in c if j - i >= 6)))
    X = fold(L, c); r = report(X, c); unpred(r, c); segdist(r)
    np.save(tag + '.npy', r)
    return c, r
