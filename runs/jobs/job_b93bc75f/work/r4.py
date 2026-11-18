from lib import *
import pickle, sys, json
L = 89
tier = {}
def add(pairs, t):
    for a, b in pairs:
        a, b = min(a, b), max(a, b)
        if b - a >= 3 and (a, b) not in tier: tier[(a, b)] = t
        elif b - a >= 3: tier[(a, b)] = min(tier[(a, b)], t)

def ladder(S, lo, hi, t0):
    """antiparallel pairs i+j=S for i in lo..hi; direct pairs tier t0, i+j=S+-2 tier t0+1, S+-1 tier t0+2"""
    for i in range(lo, hi + 1):
        add([(i, S - i)], t0)
        if i + 1 <= hi: add([(i, S - i - 2 + 0)], t0 + 1) if False else None
    for i in range(lo, hi + 1):
        j = S - i
        if i + 2 <= hi: add([(i, j), (i + 2, j), (i, j - 2)], t0 + 1)
        if i + 1 <= hi: add([(i + 1, j), (i, j - 1)], t0 + 2)

# strands
ladder(85, 33, 40, 1)     # b1-b2, turn 42-43
ladder(135, 60, 66, 1)    # b3-b4, turn 67-68
ladder(111, 45, 54, 2)    # b2-b3 long hairpin, turn 55-56
add([(41, 44), (31, 54), (32, 53), (31, 53), (31, 52), (32, 52), (59, 76), (59, 75), (60, 76), (44, 67), (45, 67)], 4)
# helices and local
add(helix(4, 12) + helix(14, 26) + helix(78, 89), 3)
add([(24, 28), (25, 28), (26, 29), (27, 30), (26, 30), (1, 4), (2, 5), (12, 15), (11, 14), (13, 16), (13, 17), (75, 78), (76, 79), (77, 80), (77, 81)], 5)
# helix 2 on face A of b1/b2
add([(26, 31), (23, 33), (22, 33), (23, 35), (22, 35), (19, 35), (19, 48), (18, 37), (25, 31), (25, 33), (26, 33)], 5)
add([(23, 50), (18, 35), (15, 37), (15, 46), (19, 50), (23, 31), (22, 48), (26, 50), (19, 46), (20, 35), (16, 37),
     (23, 48), (26, 52), (26, 54), (24, 35), (19, 37), (30, 54)], 6)
# C-term helix on face B of b2/b3/b4
add([(79, 74), (79, 73), (82, 71), (83, 71), (83, 64), (86, 69), (86, 66), (83, 73), (82, 73), (86, 68), (87, 68), (80, 74)], 5)
add([(83, 62), (86, 64), (87, 69), (86, 47), (83, 47), (89, 66), (89, 45), (79, 62), (80, 73),
     (78, 60), (79, 60), (88, 68), (85, 69), (86, 45), (82, 62), (90 - 1, 68), (87, 66), (79, 71)], 6)
# helix 1 against helix 2 and face A of b3/b4
add([(6, 19), (7, 23), (12, 18), (10, 19), (12, 19), (6, 22), (7, 22), (7, 26), (3, 26), (10, 22)], 6)
add([(6, 61), (7, 61), (9, 63), (10, 72), (6, 72), (10, 63), (3, 61), (6, 63), (9, 72)], 7)
c = sorted(tier)
print(len(c), 'long', sum(1 for a, b in c if b - a >= 6))
X = fold(L, c); r = report(X, c)
print(missed(r, c, top=60))
segs = dict(h1=(4, 12), h2=(14, 26), b1=(31, 38), b2=(44, 51), b3=(60, 66), b4=(69, 76), h3=(78, 88))
segdist(r, segs)
pickle.dump(tier, open('t4.pkl', 'wb')); np.save('r4.npy', r)
