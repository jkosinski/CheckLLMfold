from lib import *
import json
c = [tuple(x) for x in json.load(open('c3.json'))]
# round 4: drop pairs stuck at the limit that compete with better-supported packing; add a few consistently-close pairs
drop = {(4, 85), (7, 52), (66, 95), (29, 107), (77, 88), (55, 86), (22, 107), (64, 96), (2, 45), (10, 40), (13, 35), (15, 32)}
add = [(27, 103), (31, 98), (24, 107), (2, 44), (4, 52), (70, 81), (71, 82), (77, 84), (35, 90), (21, 103)]
c4 = clean([p for p in c if p not in drop] + add); print(len(c4))
X = fold(L, c4); r = report(X, c4)
json.dump(c4, open('c4.json', 'w'))
np.save('r4.npy', r)

# ranking of |i-j| >= 6 pairs
def seg(i):
    for name, a, b in [('H1', 3, 23), ('H2', 24, 45), ('B1', 50, 58), ('B2', 65, 73), ('H3', 74, 82), ('H4', 83, 108)]:
        if a <= i <= b: return name
    return 'L'
base = {('B1', 'B2'): 5.0, ('H1', 'H2'): 4.0, ('H2', 'H4'): 3.2, ('H1', 'H4'): 2.8, ('B2', 'H3'): 3.0, ('B1', 'H4'): 2.5,
        ('B2', 'H4'): 2.3, ('H2', 'B1'): 2.6, ('H3', 'H4'): 2.6, ('H1', 'B1'): 1.8}
hyd = set('AVILMFWYC')
core = {(55, 69), (53, 71), (57, 67), (51, 73), (54, 70), (56, 68), (52, 72), (21, 27), (18, 30), (14, 34), (11, 37), (7, 41), (17, 30), (10, 37), (14, 33)}
sc = []
for i, j in c4:
    if j - i < 6: continue
    s = base.get((seg(i), seg(j)), 1.5)
    s += 0.4 * (seq[i - 1] in hyd) + 0.4 * (seq[j - 1] in hyd)
    s += 0.8 * ((i, j) in core)
    s -= 0.25 * max(0.0, r[i - 1, j - 1] - 6.0)
    sc.append((-s, i, j))
sc.sort()
pairs = [(i, j) for _, i, j in sc]
print(len(pairs))
out = {"rounds": 4, "pairs": [{"i": i, "j": j} for i, j in pairs[:406]], "final_contacts": [list(p) for p in c4]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(len(out['pairs']))
