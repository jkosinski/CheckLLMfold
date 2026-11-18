from lib import *
import json, sys
L = 89
c = []
c += helix(5, 12)
c += helix(14, 27)
c += helix(75, 89)
c += [(12, 15), (11, 14), (27, 30), (28, 31), (40, 43), (41, 44), (51, 54), (66, 69), (72, 75), (73, 76), (1, 4), (2, 5)]
# hairpin b1 31-39 / b2 44-52
c += ladder_anti(31, 52, 9)
c += [(33, 48), (35, 46), (34, 51), (36, 49), (38, 47)]
# hairpin b2 45-51 / b3 54-60
c += ladder_anti(45, 60, 7)
c += [(46, 57), (48, 55), (47, 60), (49, 58)]
# hairpin b4 61-66 / b5 69-74
c += ladder_anti(61, 74, 6)
c += [(63, 70), (64, 73), (62, 71)]
# helix 14-26 on face B of b1/b2 (Y34 V36 A38 / V45 F47 V51), antiparallel to b1
c += [(26, 32), (26, 34), (26, 51), (23, 34), (23, 49), (22, 34), (23, 51), (19, 36), (19, 47), (18, 36), (15, 38), (15, 45), (19, 45), (18, 38), (22, 36), (19, 49)]
# C-term helix on face B over b3/b2 (F54 W58 Y60 / F47 V51), running back along b3
c += [(76, 60), (76, 62), (79, 58), (80, 58), (79, 60), (83, 56), (83, 54), (83, 58), (83, 47), (86, 54), (87, 54), (86, 51), (80, 47), (87, 51)]
# helix 1 / helix C side by side
c += [(17, 78), (21, 78), (21, 82), (24, 82), (24, 85), (25, 85), (24, 89), (20, 79), (17, 75), (28, 89)]
# N-term helix over the face-B residues of b4/b5 (V64 F66 I69 A71) and start of helix C
c += [(12, 60), (12, 64), (9, 64), (8, 66), (6, 66), (7, 69), (6, 71), (9, 71), (10, 76), (6, 75), (13, 62), (9, 73)]
c = clean(c, L)
print(len(c))
X = fold(L, c); r = report(X, c)
print(near(r, c)[:80])
np.save(sys.argv[0].replace('.py', '_r.npy'), r)
json.dump(c, open(sys.argv[0].replace('.py', '_c.json'), 'w'))
