from lib import *
import pickle

c = []
# helices
c += helix(25, 45)          # helix A
c += helix(76, 82)          # helix B
c += helix(84, 106)         # helix C
# hairpin strand 1 (5-10) / strand 2 (13-18), turn at 11-12
c += anti(5, 18, 6)
# strand 2 (13-17) parallel with strand 3 (51-55)
c += par(13, 51, 5)
# strand 3 (52-57) antiparallel with strand 4 (68-73)
c += anti(52, 73, 6)
# helix A over one face of strands 2/3
c += [(26, 57), (29, 55), (30, 55), (30, 57), (33, 53), (33, 55), (34, 53), (37, 53), (37, 51), (38, 51), (41, 51),
      (29, 17), (30, 17), (32, 15), (33, 15), (33, 17), (36, 15), (36, 13), (37, 13), (26, 19), (29, 19)]
# helix C over the other face of strands 3/4 (+ strand 2)
c += [(85, 73), (88, 73), (89, 73), (88, 52), (92, 54), (92, 71), (91, 52), (95, 54), (96, 56), (96, 69),
      (99, 56), (99, 69), (100, 58), (103, 58), (95, 14), (99, 14), (99, 16), (103, 16)]
# helix B against end of strand 4 and start of helix C
c += [(73, 78), (73, 77), (72, 77), (78, 85), (81, 85), (78, 88), (81, 88)]
# loops
c += [(18, 21), (19, 22), (21, 25), (22, 25), (45, 48), (46, 49), (47, 51), (58, 61), (60, 63), (62, 65), (64, 67),
      (58, 66), (59, 66), (73, 76), (80, 83), (81, 84), (2, 6), (3, 7), (104, 107), (105, 108)]
c = clean(c)
print('n contacts', len(c), 'target', 2 * L, '-', 3 * L)
X = fold(L, c)
r = report(X, c)
segs = dict(b1=(5, 10), b2=(13, 18), hA=(25, 45), b3=(51, 58), b4=(68, 73), hB=(76, 82), hC=(84, 106))
seg_dist(r, segs)
print('close but not predicted:', near_missing(r, c))
pickle.dump((c, r), open('r1.pkl', 'wb'))
