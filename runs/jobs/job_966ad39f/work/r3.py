from lib import *
import pickle, sys

def contacts3():
    c = []
    c += helix(25, 45)          # helix A (G31 in the middle, kept continuous)
    c += helix(76, 82)          # helix B
    c += helix(84, 106)         # helix C
    c += [(79, 83), (80, 83), (80, 84), (81, 84), (82, 85), (82, 86), (81, 85)]   # B -> C junction around P83
    c += anti(4, 19, 7)         # strand 1 / strand 2 hairpin, extended by one pair (model had 5-19 at 4 A)
    c += [(9, 12), (10, 14), (11, 14)]
    c += par(13, 51, 5)         # strand 2 / strand 3 parallel
    c += [(12, 50), (12, 51)]
    c += anti(52, 73, 7)        # strand 3 / strand 4 antiparallel, extended to 58-67 (model had it at 4 A)
    c += [(58, 61), (59, 62), (60, 63), (61, 64), (62, 65), (64, 67), (59, 66), (58, 66), (57, 66)]
    # helix A -> sheet, now using its apolar positions only: 26/29/30 in first half, 32/35/38/39/42 in second half
    c += [(26, 57), (29, 55), (30, 55), (30, 57), (29, 57), (32, 55), (32, 53), (35, 53), (36, 53), (35, 51), (38, 51), (39, 51), (42, 51), (42, 49),
          (29, 17), (28, 17), (32, 17), (32, 15), (35, 15), (36, 15), (35, 13), (39, 13), (38, 13), (26, 19), (29, 19), (25, 19)]
    # helix A -> strand 1 (kept only the aromatic/apolar ones the model supported)
    c += [(6, 28), (6, 29), (6, 32), (8, 32), (8, 35), (10, 35), (10, 39)]
    # helix C over the other face of strands 3/4 (+ strand 2)
    c += [(85, 73), (88, 73), (89, 73), (88, 52), (92, 54), (92, 71), (91, 52), (95, 54), (96, 56), (96, 69),
          (99, 56), (99, 69), (100, 58), (103, 58), (95, 14), (99, 14), (99, 16), (103, 16), (72, 92), (73, 91), (73, 92)]
    # helix B against the end of strand 4 and the first turns of helix C
    c += [(73, 76), (73, 77), (73, 78), (72, 77), (74, 77), (75, 78), (74, 78), (71, 78), (69, 81), (71, 81),
          (78, 85), (77, 88), (78, 88), (81, 88), (81, 89), (74, 85), (74, 88)]
    # N-terminal tail
    c += [(1, 28), (2, 28), (3, 6), (2, 6), (3, 7), (3, 20)]
    # loop 59-67 apolar residues: L66 toward the end of helix C, F64 toward start of helix A
    c += [(66, 100), (66, 103), (66, 104), (64, 26), (64, 57)]
    # loop 18-23 and the last turn of helix C; helix A N-cap; loop 46-50; C-terminus
    c += [(18, 103), (19, 104), (18, 99), (18, 21), (19, 22), (20, 23), (21, 25), (22, 25), (22, 26), (23, 26),
          (45, 48), (46, 49), (47, 50), (47, 51), (48, 51), (44, 49), (104, 107), (105, 108), (103, 107)]
    return clean(c)

if __name__ == '__main__':
    c = contacts3()
    print('n contacts', len(c), 'target', 2 * L, '-', 3 * L)
    X = fold(L, c)
    r = report(X, c)
    segs = dict(Nt=(1, 4), b1=(5, 10), b2=(13, 18), hA=(25, 45), b3=(51, 58), lp=(59, 67), b4=(68, 73), hB=(76, 82), hC=(84, 106))
    seg_dist(r, segs)
    print('close but not predicted:', near_missing(r, c, n=60))
    S = set(c)
    far = sorted([(i, j, round(float(r[i-1, j-1]), 1)) for i, j in c if r[i-1, j-1] > 7.5 and j - i >= 6], key=lambda t: -t[2])
    print('predicted long-range contacts the model holds only at the limit:', far)
    pickle.dump((c, r), open('r3.pkl', 'wb'))
