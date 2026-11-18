from lib import *
import pickle
L = 89
c = []
c += helix(4, 12) + helix(14, 26) + helix(78, 88)
c += anti(31, 52, 9)     # b1-b2 hairpin, turn 40-43
c += anti(45, 66, 7)     # b2-b3
c += anti(60, 75, 7)     # b3-b4, turn 67-68
# helix 14-26 packing on the strands
c += [(18, 33), (18, 35), (22, 33), (21, 35), (25, 31), (24, 48), (21, 46), (17, 37)]
# C-term helix on strands
c += [(83, 72), (83, 62), (86, 64), (79, 74), (82, 72), (86, 62)]
c += [(6, 19), (7, 23), (12, 18)]
c = uniq(c)
print(len(c))
X = fold(L, c); r = report(X, c)
print(missed(r, c))
segs = dict(h1=(4, 12), h2=(14, 26), b1=(31, 38), b2=(44, 51), b3=(60, 66), b4=(69, 76), h3=(78, 88))
segdist(r, segs)
pickle.dump(c, open('c1.pkl', 'wb'))
