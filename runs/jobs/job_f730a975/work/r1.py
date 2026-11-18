from lib import *
import pickle
c = []
c += anti(39, 8, 17)            # strand 8-17 with 22-31, turn at 18-21
c += helix(32, 53)
c += par(34, 22, 30)            # 22-30 with 56-64, parallel
c += anti(142, 59, 66)          # 59-66 with 76-83
c += anti(174, 77, 83)          # 77-83 with 91-97
# helix 32-53 lying across strands 22-31 and 56-65
for h, s in [(35,29),(35,30),(38,28),(39,27),(42,26),(42,27),(46,24),(46,25),(49,23),(50,22),(53,22)]:
    c.append((s, h))
for h, s in [(35,64),(35,65),(39,62),(39,63),(42,60),(42,61),(46,58),(46,59),(50,57),(50,58),(53,56)]:
    c.append((h, s))
c = clean(c); print(len(c))
X = fold(L, c); r = report(X, c)
print(missed(r, c))
segs = dict(b1=(6,17), b2=(22,31), h=(32,53), b3=(56,66), b4=(76,83), b5=(91,97), ct=(98,107), nt=(1,5))
segdist(r, segs)
pickle.dump(c, open('c1.pkl','wb')); np.save('X1.npy', X)
