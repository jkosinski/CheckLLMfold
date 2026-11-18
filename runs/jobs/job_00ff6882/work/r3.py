from lib import *
import json
c = [tuple(x) for x in json.load(open('c2.json'))]
# round 3: drop speculative strand1-helix1 pairs the model could not fit together with the helix-helix packing
drop = {(14, 56), (10, 54), (6, 44), (58, 95), (64, 92), (78, 88)}
# add: C-terminal halves of helix 1 / helix 4 and N-term of helix 2 were under-packed (low contact number)
add = [(15, 98), (18, 98), (19, 102), (22, 102), (22, 106), (19, 103), (26, 102), (26, 106), (29, 103), (30, 98), (27, 98),
       (22, 107), (25, 107), (26, 107), (29, 107), (23, 106), (15, 95), (18, 102), (30, 94), (34, 91), (15, 29), (19, 28), (21, 28), (15, 32), (11, 36)]
# N-terminal tail / strand 1 / helix 2 end pairs seen close in model and consistent with hydrophobic pattern (Y6, V7 with F54, V44, L38)
add += [(7, 54), (6, 54), (4, 53), (5, 41), (8, 40), (7, 37), (41, 53), (35, 56), (54, 88), (56, 91), (54, 87), (71, 84), (72, 85), (77, 86), (56, 66)]
c3 = clean([p for p in c if p not in drop] + add); print(len(c3))
X = fold(L, c3); r = report(X, c3)
print(missed(r, c3, thr=5.5))
json.dump(c3, open('c3.json', 'w'))
cn = (r < 8).sum(1) - 1
print(''.join(str(min(9, x // 3)) for x in cn)); print(seq)
