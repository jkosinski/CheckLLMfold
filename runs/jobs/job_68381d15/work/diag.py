import sys, json
from lib import *
name = sys.argv[1]
cl = json.load(open(name + '.contacts.json'))
contacts = [(a, b) for a, b, t in cl]
L = 107
X = fold(L, contacts)
c = X.mean(0)
d = np.linalg.norm(X - c, axis=1)
print('distance from centroid per residue:')
for s in range(0, L, 10):
    print(s + 1, ' '.join('%4.0f' % v for v in d[s:s + 10]))
r = np.linalg.norm(X[:, None] - X[None], axis=-1)
nb = (r < 10).sum(1) - 1
print('neighbours within 10 A:')
for s in range(0, L, 10):
    print(s + 1, ' '.join('%4d' % v for v in nb[s:s + 10]))
segs = {'N': (1, 5), 'A': (7, 17), 'B': (22, 31), 'H': (33, 52), 'C': (57, 66), 'L1': (67, 75), 'D': (76, 82), 'L2': (83, 89), 'E': (90, 97), 'T': (98, 107)}
names = list(segs)
print('segment centroid distances:')
cen = {k: X[a - 1:b].mean(0) for k, (a, b) in segs.items()}
print('     ' + ' '.join('%5s' % n for n in names))
for a in names:
    print('%5s' % a + ' ' + ' '.join('%5.1f' % np.linalg.norm(cen[a] - cen[b]) for b in names))
