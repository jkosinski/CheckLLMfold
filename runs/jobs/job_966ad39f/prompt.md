We are testing whether you can tell, from a protein sequence alone, which residues are in contact in its folded structure. You will predict a contact map, turn it into a 3D model with the code below, look at what the model says about your contacts, and revise them, several times.

Rules
- Predict contacts from the sequence. You may reason about secondary structure (helices, strands, turns) from the sequence and from your local contacts. Do not reason about domains, the protein's family or what it is, or its overall fold or topology, and do not describe any of these. Do not search the web or any database. What we want to measure is whether you can see contacts in the sequence itself.
- Residues are numbered from 1. A contact means C-beta atoms (C-alpha for glycine) within 8 A, at sequence separation |i - j| >= 3.

Sequence

MNHEDYVQSVHGQIVKTAMAMQSGEMSYLLGARKLDALRHEASVKNDDEDFMVFVAIASDTDDFPLGTVREYWDKDALAKLQPEIDAAETWAREHSEHICAKLIRRFS

Procedure

1. Predict contacts. Write a list of (i, j) pairs you believe are in contact, local ones (|i - j| from 3 to 5) included. A folded protein has about 3 x L of them; give between 2 x L and 3 x L.
2. Build a model. Run the code below once in your code tool, then:
   X = fold(L, contacts)
   r = report(X, contacts)
   fold() places consecutive residues 3.8 A apart and every contact at about 6 A, fills in the rest as shortest paths, embeds that in 3D, and refines it so contacts stay within 8 A and no two residues overlap. Given the true contacts it rebuilds native structures at TM-score 0.64-0.84, so the model is only as good as your contacts.
3. Revise. Use the report and the model:
   - contacts it could not satisfy contradict the others: decide which are wrong;
   - a radius of gyration well above the compact value means the chain is not held together: you are missing contacts;
   - r is the model's distance matrix: pairs your model puts close together that you did not predict may be contacts you missed, and contacts you did predict may now look implausible.
   Rewrite the contact list and rebuild. Do at least 3 and at most 6 rounds, each time saying in a sentence or two which contacts you changed and why (in terms of contacts, secondary structure and the model, not of folds or protein families).
4. Answer. From your final contacts, rank the 406 pairs with |i - j| >= 6 you are most confident of, most confident first.
Also give your final-round contact list (all separations, |i - j| >= 3) in the JSON as "final_contacts".

Return, at the end, exactly one JSON object:
{"rounds": <number of rounds>, "pairs": [{"i": <int>, "j": <int>}, ...], "final_contacts": [[<int>, <int>], ...]}

Code (numpy only):

```python
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
    """contacts: list of (i, j), 1-based, |i - j| >= 3. Consecutive residues 3.8 A, each contact 6 A, every other
    distance the shortest path through those; embedded in 3D by classical MDS. Returns (X [L, 3], distances)."""
    d = np.full((L, L), np.inf); np.fill_diagonal(d, 0.0)
    k = np.arange(L - 1); d[k, k + 1] = d[k + 1, k] = 3.8
    for i, j in contacts: d[i - 1, j - 1] = d[j - 1, i - 1] = 6.0
    return classic_mds(floyd_warshall(d)), d

def refine(X, contacts, steps=3000, lr=0.01):
    """gradient descent: consecutive residues 3.8 A apart, every contact within 8 A, no two residues closer than 4 A"""
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
    """MDS cannot tell a fold from its mirror image. Where four consecutive residues are wound tightly (CA i..i+3
    under 6 A, as in a helix) the turn should be right-handed (dihedral ~ +50 deg); mirror the model if most are not."""
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
```
