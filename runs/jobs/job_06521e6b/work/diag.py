import sys, numpy as np
from lib import *
X = np.load(sys.argv[1]); c = X.mean(0); d = np.linalg.norm(X - c, axis=1)
r = np.linalg.norm(X[:, None] - X[None], axis=-1)
nb = (r < 10).sum(1) - 1
for a in range(0, L, 12):
    print(' '.join('%s%d:%.0f/%d' % (SEQ[i], i + 1, d[i], nb[i]) for i in range(a, min(L, a + 12))))
segs = {'b1': (5, 10), 'a1': (13, 22), 'a2': (25, 42), 'b2': (51, 58), 'b3': (66, 73), 'a3': (75, 82),
        'a4a': (84, 95), 'a4b': (96, 106)}
cen = {k: X[a - 1:b].mean(0) for k, (a, b) in segs.items()}
ks = list(segs)
print('     ' + ' '.join('%5s' % k for k in ks))
for k in ks: print('%5s' % k + ' '.join('%5.1f' % np.linalg.norm(cen[k] - cen[m]) for m in ks))
