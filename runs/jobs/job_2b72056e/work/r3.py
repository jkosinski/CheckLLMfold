from lib import *
from r2 import r2list

def r3list():
    c = set((min(a, b), max(a, b)) for a, b in r2list())
    # residue 30 cannot reach both strand 11-17 and strand 40-45 (two strands apart): keep the helix C-term on 5/69/43/41
    # E16 points away from the helix face; 17-27, 1-63, 50-56 were stretched to the limit
    for p in [(13, 30), (15, 30), (16, 23), (17, 27), (1, 63), (50, 56), (20, 23)]:
        c.discard(p)
    # pairs the model brings close that fit the secondary structure
    c |= {(34, 40), (37, 72), (37, 73), (36, 72), (38, 73),      # loop 35-39 running beside the end of strand 66-72
          (45, 60), (46, 60), (60, 66), (49, 59), (47, 59), (44, 59),  # 56-61 turn docking on the 44-49 hairpin
          (1, 22), (2, 22), (17, 24), (15, 25),                  # loop 18-22 / helix N-term over strand ends 1-3, 15-17
          (27, 69), (26, 30 + 13), (34, 38), (31, 36)}
    return sorted(c)

if __name__ == '__main__':
    contacts, X, r = run(r3list(), tag='R3')
    np.save('r3.npy', r)
