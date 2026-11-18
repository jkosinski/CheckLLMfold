from lib import *
from ideal import *
import json, sys
L = 89
P = layout()
P2 = {}
for k, (x, y, z) in P.items():
    if k <= 27 or k >= 70: z = z - 1.0 if z > 0 else z + 1.0
    P2[k] = (x, y, z)
P = P2
CUT = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
ideal = contacts_from(P, cut=CUT)
allid = contacts_from(P, cut=30.0, sep=6)
c = []
c += helix(5, 12) + helix(14, 27) + helix(70, 89)
c += [(1, 4), (2, 5), (11, 14), (12, 15), (27, 30), (28, 31), (26, 31), (40, 43), (41, 44), (40, 44), (51, 54), (50, 54), (50, 55), (51, 55),
      (58, 61), (57, 60), (58, 62), (57, 61), (66, 69), (67, 70), (66, 70), (65, 69)]
c += list(ideal)
c += [(19, 35), (23, 33), (6, 56), (4, 21), (3, 24), (26, 33), (25, 31), (27, 32),
      (30, 52), (31, 53), (69, 72), (69, 73), (66, 73), (64, 72), (62, 76), (60, 79), (60, 83), (58, 79),
      (45, 87), (38, 89), (36, 86), (38, 86), (45, 83), (53, 72), (51, 79)]
c = clean(c, L)
nl = sum(1 for a, b in c if b - a >= 6)
print(len(c), 'long (>=6):', nl)
X = fold(L, c); r = report(X, c)
print(near(r, c)[:40])

# ---------- ranking of pairs with |i-j| >= 6 ----------
seg = lambda a: 'N' if a <= 4 else 'hN' if a <= 12 else 'h1' if a <= 27 else 'l' if a <= 30 else 'b1' if a <= 40 else 'b2' if a <= 53 else 'b3' if a <= 59 else 'b4' if a <= 66 else 'hC'
S = set(c)
def conf(a, b):
    sa, sb = seg(a), seg(b); key = sa + '-' + sb
    di = allid.get((a, b)); dm = float(r[a - 1, b - 1])
    d = di if di is not None else dm + 1.0
    base = {'b1-b2': 1.00, 'b2-b3': 0.93, 'b3-b4': 0.86, 'h1-b1': 0.62, 'h1-b2': 0.52, 'l-b2': 0.60, 'h1-l': 0.5,
            'b4-hC': 0.60, 'b3-hC': 0.55, 'b2-hC': 0.48, 'b1-hC': 0.42, 'hN-h1': 0.45, 'N-h1': 0.35,
            'b1-b3': 0.40, 'b2-b4': 0.40}.get(key, 0.30)
    s = base - 0.035 * d
    if sa == sb: s -= 0.3            # i,i+6.. inside one helix/strand is not a contact
    if (a, b) in S: s += 1.0         # final contacts first
    return s
pairs = [(a, b) for a in range(1, L + 1) for b in range(a + 6, L + 1)]
pairs.sort(key=lambda p: -conf(*p))
top = pairs[:335]
print('final contacts among top335:', sum(1 for p in top if p in S), 'of', nl)
print(top[:40]); print(top[nl - 5:nl + 30]); print(top[-20:])
out = {'rounds': 6, 'pairs': [{'i': a, 'j': b} for a, b in top], 'final_contacts': [[a, b] for a, b in c]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
