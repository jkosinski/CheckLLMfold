from lib import *
from ideal import *
import json, sys
L = 89
P = layout()
ideal = contacts_from(P, cut=8.0)
c = []
c += helix(5, 12) + helix(14, 27) + helix(70, 89)
c += [(1, 4), (2, 5), (11, 14), (12, 15), (27, 30), (28, 31), (26, 31), (40, 43), (41, 44), (40, 44), (51, 54), (58, 61), (57, 60), (66, 69), (67, 70), (66, 70), (65, 69)]
# sheet + helix packing from the idealised layout (helix-internal pairs come from helix())
for (a, b), d in ideal.items():
    c.append((a, b))
c = clean(c, L)
print(len(c), 'long (>=6):', sum(1 for a, b in c if b - a >= 6))
X = fold(L, c); r = report(X, c)
print(near(r, c)[:70])
np.save(sys.argv[0].replace('.py', '_r.npy'), r)
json.dump(c, open(sys.argv[0].replace('.py', '_c.json'), 'w'))
seg = lambda a: 'hN' if a <= 12 else 'h1' if a <= 27 else 'b1' if a <= 40 else 'b2' if a <= 52 else 'b3' if a <= 58 else 'b4' if a <= 66 else 'hC'
from collections import defaultdict
g = defaultdict(list)
for a, b in c:
    if b - a >= 6 and seg(a) != seg(b): g[seg(a) + '-' + seg(b)].append((a, b))
for k, v in g.items(): print(k, len(v), v)
