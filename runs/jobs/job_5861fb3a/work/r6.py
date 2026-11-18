from lib import *
from ideal import *
import json, sys
L = 89
P = layout()
# pull the helices 1 A closer to the sheet than in round 5
P2 = {}
for k, (x, y, z) in P.items():
    if k <= 27 or k >= 70: z = z - 1.0 if z > 0 else z + 1.0
    P2[k] = (x, y, z)
P = P2
ideal = contacts_from(P, cut=8.6)
allid = contacts_from(P, cut=14.0, sep=6)
c = []
c += helix(5, 12) + helix(14, 27) + helix(70, 89)
c += [(1, 4), (2, 5), (11, 14), (12, 15), (27, 30), (28, 31), (26, 31), (40, 43), (41, 44), (40, 44), (51, 54), (50, 54), (50, 55), (51, 55),
      (58, 61), (57, 60), (58, 62), (57, 61), (66, 69), (67, 70), (66, 70), (65, 69)]
c += list(ideal)
# pairs the round-5 model put close that fit the same packing; N-terminal tail next to the 14-26 helix; helix-1 end to strand start
c += [(19, 35), (23, 33), (10, 60), (6, 56), (4, 21), (2, 8), (3, 24), (3, 21), (1, 24), (2, 25), (26, 33), (25, 31), (27, 32), (23, 31),
      (30, 52), (30, 53), (31, 53), (32, 53), (69, 72), (69, 73), (66, 73), (64, 72), (62, 76), (60, 79), (60, 83), (58, 79), (58, 80),
      (45, 87), (38, 89), (37, 89), (45, 89), (36, 86), (38, 86), (45, 83), (53, 72), (53, 75), (51, 79), (52, 75)]
c = clean(c, L)
print(len(c), 'long (>=6):', sum(1 for a, b in c if b - a >= 6))
X = fold(L, c); r = report(X, c)
print(near(r, c)[:70])
np.save('r6_r.npy', r)
json.dump(c, open('r6_c.json', 'w'))
json.dump({'%d,%d' % k: v for k, v in allid.items()}, open('r6_ideal.json', 'w'))
