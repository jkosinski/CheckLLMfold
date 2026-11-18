from lib import *
import pickle, json
import numpy as np
cs = set(pickle.load(open('c3.pkl', 'rb')))
# round 4: drop borderline contacts the model keeps stretching to exactly 8 A and that rest on weak reasoning;
# add the helix-start / strand-start and 45-48 turn / 60-62 contacts the model places close.
for p in [(13, 30), (23, 55), (36, 71), (45, 59)]: cs.discard(p)
cs |= {(1, 23), (15, 25), (45, 62), (48, 62), (46, 61), (49, 61), (1, 61)}
cs, X, r = run(cs, 'round4'); pickle.dump(cs, open('c4.pkl', 'wb'))
for p in [(13, 30), (23, 55), (36, 71), (45, 59), (17, 61)]: print(p, round(float(r[p[0]-1, p[1]-1]), 1))

# ---- confidence ranking of pairs with |i-j| >= 6 ----
score = {}
def add(pairs, s):
    for i, j in norm(pairs):
        if j - i >= 6: score[(i, j)] = max(score.get((i, j), 0), s)
add(anti(18, 1, 7), 10.0)                      # hairpin 1-7 / 11-17 register
add(anti(112, 40, 46), 9.5)                    # 40-46 / 66-72 register
add(par(64, 1, 7), 9.0)                        # 1-7 / 65-71 parallel register
add(anti(18, 1, 7, (-1, 1)), 8.5)
add(anti(112, 40, 46, (-1, 1)), 8.0)
add(par(64, 1, 7, (-1, 1)), 7.5)
add(anti(18, 1, 7, (-2, 2)), 7.0)
add(anti(93, 42, 45, (-1, 0, 1)), 7.0)
add([(42, 72), (44, 70), (41, 69), (43, 67), (2, 64), (1, 63)], 6.5)
helix_sheet = [(3, 23), (17, 23), (15, 23), (3, 26), (5, 26), (15, 26), (5, 30), (15, 30), (26, 43), (30, 43), (30, 41),
               (27, 43), (34, 41), (26, 67), (30, 69), (26, 69), (5, 67), (5, 69), (3, 67), (16, 23), (1, 23)]
add(helix_sheet, 6.0)
add([(23, 56), (26, 56), (61, 67), (3, 61), (45, 61), (50, 59), (50, 56), (43, 56), (56, 67), (27, 56), (22, 56)], 5.5)
final = set(cs)
ranked = []
for (i, j) in final:
    if j - i >= 6:
        s = score.get((i, j), 4.5)
        ranked.append((-(s - 0.05 * float(r[i-1, j-1])), i, j))
ranked.sort()
pairs = [(i, j) for _, i, j in ranked]
# fill to 286 with the next-closest pairs in the final model
rest = sorted([(float(r[i-1, j-1]), i, j) for i in range(1, 77) for j in range(i + 6, 77) if (i, j) not in final])
for d, i, j in rest:
    if len(pairs) >= 286: break
    pairs.append((i, j))
print('n final contacts', len(cs), 'n long predicted', len(ranked), 'n ranked', len(pairs), 'last fill distance', round(d, 1))
assert len(pairs) == 286 and len(set(pairs)) == 286
out = {"rounds": 4, "pairs": [{"i": i, "j": j} for i, j in pairs], "final_contacts": [[i, j] for i, j in cs]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(json.dumps(out, separators=(',', ':')))
