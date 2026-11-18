from lib import *
import pickle, sys
c = []
c += anti(39, 6, 17)            # strand 6-17 with 22-33 (extended by two pairs at the open end)
c += helix(32, 53)
c += par(34, 22, 31)            # 22-31 with 56-65, parallel
c += anti(142, 59, 66)          # 59-66 with 76-83
c += anti(174, 77, 83)          # 77-83 with 91-97
# helix on the side of 23/25/27/29 and 57/59/61/63/65
hs = [(35,29),(35,27),(38,29),(39,27),(42,27),(42,25),(46,25),(46,23),(49,23),(50,23),
      (35,63),(35,65),(38,63),(39,63),(39,61),(42,61),(42,59),(46,59),(46,57),(49,57),(50,57),(53,57),(53,56)]
c += hs
# loop 67-75 closes on itself, 84-90 turn
c += [(67,75),(68,74),(69,73),(66,75),(84,90),(85,89),(83,90)]
# tails: 1-5 alongside the helix start; 98-107 folds back on 91-97 and 77-83
c += [(4,34),(5,34),(3,36),(4,37),(2,37)]
c += [(97,100),(96,101),(95,102),(94,103),(93,104),(92,105),(96,100),(94,102),(92,104),(81,104),(79,106),(83,102),(80,105),(78,107)]
c = clean(c); print(len(c))
X = fold(L, c); r = report(X, c)
print(missed(r, c))
segs = dict(b1=(6,17), b2=(22,31), h=(32,53), b3=(56,66), b4=(76,83), b5=(91,97), ct=(98,107), nt=(1,5), lp=(67,75))
segdist(r, segs)
pickle.dump(c, open('c2.pkl','wb')); np.save('X2.npy', X)
