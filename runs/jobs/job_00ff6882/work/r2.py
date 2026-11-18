from lib import *
import json, sys
c = [tuple(x) for x in json.load(open('c1.json'))]
X = fold(L, c); r = np.linalg.norm(X[:, None] - X[None], axis=-1)
# round 2: accept model-close unpredicted pairs (< 6.0 A) between elements already packed,
add = [(i, j) for i, j, d in missed(r, c, thr=6.0)]
# tie N-terminal segment (1-6) to end of helix 2 / start of strand 1, and helix 3 to helix 4 start + strand 2
add += [(3, 44), (4, 41), (3, 41), (2, 45), (4, 45), (6, 44), (3, 51), (6, 53), (7, 53), (10, 55), (7, 85), (4, 85)]
add += [(77, 88), (78, 85), (78, 88), (81, 85), (74, 85), (77, 73), (80, 67), (77, 69)]
# strand 1 other face toward helix 1 / helix 4 middle
add += [(54, 10), (56, 14), (52, 7), (56, 92), (58, 95), (64, 92), (64, 96), (66, 95)]
# drop weakly supported: those sitting exactly at the limit that compete
drop = {(53, 42), (67, 81)}
c2 = clean([p for p in c if p not in drop] + add); print(len(c2))
X = fold(L, c2); r = report(X, c2)
print(missed(r, c2, thr=6.0))
json.dump(c2, open('c2.json', 'w'))
cn = (r < 8).sum(1) - 1
print(''.join(str(min(9, x // 3)) for x in cn)); print(seq)
