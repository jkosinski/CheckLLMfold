from lib import *
import json
T = {}
def add(pairs, tier, store=None):
    store = T if store is None else store
    for i, j in pairs:
        i, j = min(i, j), max(i, j)
        if j - i < 3 or not (1 <= i and j <= L): continue
        store[(i, j)] = min(tier, store.get((i, j), 99))

LADDERS = [(82, 31, 38, 1), (156, 67, 74, 1), (28, 6, 12, 1), (104, 31, 38, 2), (77, 16, 22, 2)]
for S, lo, hi, t in LADDERS:
    add(anti(S, lo, hi, (0,)), t)
    add(anti(S, lo, hi, (1, -1)), t + 2)
add(par(78, 7, 11, (0,)), 5)
add(par(78, 7, 11, (1, -1)), 6)
add([(1,4),(2,5),(3,6),(4,7),
     (12,15),(13,16),(12,16),(11,16),
     (22,25),(23,26),(24,27),(25,28),(26,29),(27,30),(24,28),(25,29),(26,30),(23,27),
     (37,40),(38,41),(39,42),(40,43),(41,44),(38,42),(39,43),(38,43),(39,44),(40,44),
     (49,52),(50,53),(51,54),(52,55),(50,54),(51,55),(49,55),
     (60,63),(61,64),(62,65),(63,66),(61,65),(62,66),(61,66),
     (73,76),(74,77),(75,78),(76,79),(77,80),(78,81),(79,82),(80,83),(75,79),(76,80),(77,81),(78,82),
     (75,82),(74,82),(74,83),(76,81),(76,82)], 4)
add([(21,33),(19,35),(21,35),(19,33),(21,71),(19,69),(21,69),(19,71),(56,33),(58,35),(58,33),(56,35),
     (56,49),(58,47),(56,47),(58,49)], 5)
add([(17,37),(21,31),(17,35),(17,67),(19,67),(60,37),(58,37),(60,45),(58,45),(7,71),(9,69),(9,71)], 6)
add([(25,33),(25,21),(28,31),(28,33),(25,49),(22,31),(25,31)], 6)
add([(37,61),(38,61),(39,61),(39,63),(40,61),(15,62),(16,66),(16,64),(13,89),(14,89),
     (15,65),(16,65),(15,66),(41,63),(39,64),(40,63)], 7)
add([(25,54),(25,55),(27,54),(28,51),(29,51),(28,53),(30,51),(30,52),(26,54),(27,53),(28,52),
     (5,23),(4,23),(5,84),(6,84),(5,73),(29,74),(30,74),(28,75),(30,73),(23,73),(22,73),(30,75)], 7)
# new in round 5
add([(60,66),(61,67),(60,67),(38,62),(39,62),(26,53),(28,50),(29,75),(12,89),(14,65),(28,32)], 7)

c = sorted(T)
c, r = run(c, 'r5')
assert 2 * L <= len(c) <= 3 * L, len(c)

# ranking: final contacts with |i-j| >= 6 by tier, then by model distance
ranked = sorted([k for k in c if k[1] - k[0] >= 6], key=lambda k: (T[k], r[k[0]-1, k[1]-1]))
print('final contacts with sep>=6:', len(ranked))
# fill up to 335: register-shifted ladder pairs, then inward-facing packing pairs, then model distance
F = {}
for S, lo, hi, t in LADDERS:
    add(anti(S, lo, hi, (2, -2)), 8, F)
inward = {'s1': [7, 9, 11], 's2': [17, 19, 21], 's5': [56, 58, 60], 's3': [31, 33, 35, 37], 's4': [45, 47, 49],
          's6': [67, 69, 71, 73], 's7': [85, 87, 89]}
for a, b in [('s2', 's3'), ('s2', 's6'), ('s5', 's3'), ('s5', 's4'), ('s1', 's6'), ('s1', 's7')]:
    add([(i, j) for i in inward[a] for j in inward[b] if r[i-1, j-1] < 7.0], 8, F)
for i in range(1, L + 1):
    for j in range(i + 6, L + 1):
        if r[i-1, j-1] < 8.0: add([(i, j)], 9, F)
fill = sorted([k for k in F if k not in T and k[1] - k[0] >= 6], key=lambda k: (F[k], r[k[0]-1, k[1]-1]))
print('fill candidates:', len(fill))
N = 335
pairs = (ranked + fill)[:N]
assert len(pairs) == N and len(set(pairs)) == N and all(j - i >= 6 for i, j in pairs)
out = {'rounds': 5, 'pairs': [{'i': i, 'j': j} for i, j in pairs], 'final_contacts': [[i, j] for i, j in c]}
s = json.dumps(out, separators=(',', ':'))
open('answer.json', 'w').write(s)
print(len(c), len(pairs))
print(s)
