from lib import *
import pickle, sys

def base():
    c = []
    c += helix(25, 45)          # helix A
    c += helix(76, 82)          # helix B
    c += helix(84, 106)         # helix C
    c += anti(5, 18, 6)         # strand 1 / strand 2 hairpin
    c += par(13, 51, 5)         # strand 2 / strand 3 parallel
    c += anti(52, 73, 6)        # strand 3 / strand 4 antiparallel
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
    return c

if __name__ == '__main__':
    c = base()
    # round 2 additions
    # helix A also reaches strand 1 (model put 6/7 next to 29-32 and 9/10 next to 34-37)
    c += [(6, 29), (6, 30), (8, 33), (8, 30), (10, 37), (10, 36), (10, 33)]
    # helix B was detached: tie it to strand 4 and to the first turns of helix C
    c += [(74, 77), (75, 78), (74, 78), (71, 78), (69, 81), (71, 81), (77, 88), (81, 89), (82, 85), (82, 86), (74, 85), (74, 88), (72, 92), (73, 91), (73, 92)]
    # N-terminal tail along the start of helix A / strand 1
    c += [(1, 28), (2, 28), (3, 6), (4, 7), (4, 19), (2, 20)]
    # loop 59-67: L66 toward the end of helix C, F64 toward start of helix A
    c += [(66, 100), (66, 103), (66, 104), (64, 26), (64, 57), (61, 64), (59, 62), (57, 66)]
    # loop 18-23 and the last turn of helix C
    c += [(18, 103), (19, 104), (21, 107), (18, 99), (20, 23), (22, 26), (23, 26)]
    # turn of hairpin 1/2 and loop 46-50
    c += [(9, 12), (10, 14), (11, 14), (48, 51), (47, 50), (12, 50), (12, 51), (44, 49), (41, 49)]
    c = clean(c)
    print('n contacts', len(c), 'target', 2 * L, '-', 3 * L)
    X = fold(L, c)
    r = report(X, c)
    segs = dict(Nt=(1, 4), b1=(5, 10), b2=(13, 18), hA=(25, 45), b3=(51, 58), lp=(59, 67), b4=(68, 73), hB=(76, 82), hC=(84, 106))
    seg_dist(r, segs)
    print('close but not predicted:', near_missing(r, c, n=60))
    pickle.dump((c, r), open('r2.pkl', 'wb'))
