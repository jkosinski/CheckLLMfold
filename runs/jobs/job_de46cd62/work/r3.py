from lib import *
import pickle, sys
import numpy as np
base = pickle.load(open('c2.pkl', 'rb'))
cs = set(base)
# residue 61 cannot sit over both strand 11-17 (17) and strand 40-45 (43): two strands lie between them.
# Keep it on the 43/45/67 side; drop the 17-side contacts of 56 and 61 and the weak 19-56.
for p in [(17, 61), (17, 56), (19, 56), (1, 61)]: cs.discard(p)
# model puts the end of strand 11-17 against the first helix turn, and 59-64 against the 45-48 turn: add those
cs |= {(16, 23), (15, 27), (16, 26), (45, 60), (46, 62), (46, 64), (44, 61), (48, 60), (49, 60), (21, 55), (21, 57),
       (26, 68), (61, 68), (1, 67), (18, 63), (50, 60), (57, 63)}
cs, X, r = run(cs, 'round3'); pickle.dump(cs, open('c3.pkl', 'wb'))
for p in [(17, 61), (17, 56), (19, 56), (1, 61)]: print(p, round(float(r[p[0]-1, p[1]-1]), 1))
np.save('r3.npy', r); np.save('X3.npy', X)
# secondary-structure sanity: CA(i)-CA(i+2), CA(i)-CA(i+3) along chain
print('i,i+3 distances:', [round(float(r[i, i+3]), 1) for i in range(73)])
