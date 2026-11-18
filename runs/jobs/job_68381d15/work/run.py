import sys, json
from lib import *

seq = "MGRAGWYGVRCVFRWVHEGRQVYEETVTVWRAESFEEALQKAEVGAFEYAAECDGQYLEFAQAYFIGEDEVIGEGAEVFSLMRESELGGRDYVTRYFDTGDERQGAW"
L = len(seq)
C = {}

def add(i, j, t):
    if i > j: i, j = j, i
    if i < 1 or j > L or j - i < 3: return
    if (i, j) not in C or C[(i, j)] > t: C[(i, j)] = t

def anti(s, lo, hi, t, width=1):
    for i in range(lo, hi + 1):
        add(i, s - i, t)
        if width >= 1:
            add(i, s - i + 1, t + 1); add(i, s - i - 1, t + 1)
        if width >= 2:
            add(i, s - i + 2, t + 2); add(i, s - i - 2, t + 2)

def para(off, lo, hi, t, width=1):
    for i in range(lo, hi + 1):
        add(i, i + off, t)
        if width >= 1:
            add(i, i + off + 1, t + 1); add(i, i + off - 1, t + 1)

def helix(a, b, t=1):
    for i in range(a, b + 1):
        for k in (3, 4):
            if i + k <= b: add(i, i + k, t)

name = sys.argv[1]
cfg = json.load(open(name))
for h in cfg.get('helix', []): helix(*h)
for a in cfg.get('anti', []): anti(*a)
for a in cfg.get('para', []): para(*a)
for i, j, t in cfg.get('pairs', []): add(i, j, t)
for i, j in cfg.get('remove', []):
    C.pop((min(i, j), max(i, j)), None)
contacts = sorted(C)
print('L', L, 'n contacts', len(contacts), 'long(>=6)', sum(1 for i, j in contacts if j - i >= 6))
X = fold(L, contacts); r = report(X, contacts)
np.save(name + '.r.npy', r)
json.dump([[i, j, C[(i, j)]] for i, j in contacts], open(name + '.contacts.json', 'w'))
P = set(contacts)
cand = [(round(float(r[i - 1, j - 1]), 1), i, j) for i in range(1, L + 1) for j in range(i + 3, L + 1)
        if (i, j) not in P and r[i - 1, j - 1] < 7]
cand.sort(); print('close (<7) but unpredicted:', len(cand)); print(cand[:90])
far = [(i, j, round(float(r[i - 1, j - 1]), 1)) for i, j in contacts if r[i - 1, j - 1] > 7.5]
print('predicted but >7.5:', far)
