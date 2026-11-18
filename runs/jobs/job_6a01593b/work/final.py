import json
import numpy as np
from lib import SEQ, fold, report
from r1 import H, core, flank
from r2 import run, norm
import r3

L = len(SEQ)
r3.ADD[5] = [(79, 90), (79, 89), (80, 86), (17, 92), (15, 94), (18, 26)]
r3.DROP[5] = [(10, 34), (3, 101), (23, 79), (72, 98)]
c = r3.make(5)
c, X, r = run(5, c)
assert 2 * L <= len(c) <= 3 * L, len(c)
assert all(1 <= i < j <= L and j - i >= 3 for i, j in c)

def helix_of(i):
    for h, (a, b) in H.items():
        if a <= i <= b:
            return h
    return 0

HAIR = {(1, 2), (3, 4)}
ADJ = {(2, 3), (1, 4)}
S = {(1, 2): 43.5, (3, 4): 166, (2, 3): 98, (1, 4): 108}

def score(i, j, predicted):
    a, b = helix_of(i), helix_of(j)
    d = float(r[i - 1, j - 1])
    if predicted:
        if a == 0 or b == 0 or a == b:
            base = 35
        elif (a, b) in HAIR or (a, b) in ADJ:
            nc = (i in core[a]) + (j in core[b])
            k = abs(i + j - S[(a, b)])
            top = 100 if (a, b) in HAIR else 80
            base = top if (nc == 2 and k <= 1) else top - 12 if k <= 1 else top - 22
            if (a, b) in HAIR:
                base += max(0, 6 - 0.25 * (j - i))  # closer to the connecting turn = more certain
        else:
            base = 42
        return base - 1.5 * d
    if a and b and (a, b) in HAIR:
        return 62 - 4 * d
    if a and b and (a, b) in ADJ:
        return 52 - 4 * d
    return 40 - 4 * d

Sc = set(c)
allp = [(score(i, j, (i, j) in Sc), i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1)]
allp.sort(reverse=True)
top = allp[:399]
print('predicted long pairs:', sum(1 for i, j in c if j - i >= 6), ' in top399:', sum(1 for _, i, j in top if (i, j) in Sc))
print('first 25:', [(i, j) for _, i, j in top[:25]])
print('last 10:', [(round(s, 1), i, j, round(float(r[i - 1, j - 1]), 1)) for s, i, j in top[-10:]])
out = {"rounds": 5, "pairs": [{"i": i, "j": j} for _, i, j in top], "final_contacts": [[i, j] for i, j in c]}
assert len(out["pairs"]) == 399 and len(set((p["i"], p["j"]) for p in out["pairs"])) == 399
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(len(json.dumps(out, separators=(',', ':'))))
