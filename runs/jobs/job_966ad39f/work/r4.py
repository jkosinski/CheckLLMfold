from lib import *
import pickle, sys, json

# contact groups, in decreasing order of confidence (used later for ranking)
def groups(round_no=4):
    g = []
    g.append(('hairpin12', [(4 + k, 19 - k) for k in range(7)]))
    g.append(('hairpin34', [(52 + k, 73 - k) for k in range(7)]))
    g.append(('par23', [(13 + k, 51 + k) for k in range(5)]))
    g.append(('hairpin12_diag', [p for p in anti(4, 19, 7) if p[0] + p[1] != 23] + [(10, 14)]))
    g.append(('hairpin34_diag', [p for p in anti(52, 73, 7) if p[0] + p[1] != 125]))
    g.append(('par23_diag', [p for p in par(13, 51, 5) if p[1] - p[0] != 38] + [(12, 50), (12, 51)]))
    hA_sheet = [(29, 17), (32, 15), (35, 15), (32, 17), (29, 55), (30, 55), (32, 53), (35, 53), (26, 57), (30, 57), (29, 57), (32, 55),
                (36, 53), (35, 51), (38, 51), (39, 51), (42, 51), (42, 49), (28, 17), (36, 15), (35, 13), (39, 13), (38, 13),
                (26, 19), (29, 19), (25, 19)]
    hC_sheet = [(92, 54), (96, 56), (99, 56), (92, 71), (96, 69), (88, 73), (95, 54), (99, 69), (85, 73), (89, 73), (88, 52), (91, 52),
                (100, 58), (103, 58), (72, 92), (73, 91), (73, 92), (95, 14), (99, 14), (99, 16), (103, 16)]
    g.append(('helixA_sheet', hA_sheet))
    g.append(('helixC_sheet', hC_sheet))
    hB = [(73, 77), (73, 78), (72, 77), (78, 85), (74, 85), (74, 88), (71, 78), (78, 88), (71, 81), (82, 86)]
    if round_no >= 5:
        hB += [(75, 85), (72, 88), (77, 85)]
    g.append(('helixB', hB))
    g.append(('helixA_strand1', [(6, 28), (6, 29), (6, 32), (8, 32), (8, 35), (10, 35), (10, 39)]))
    lp = [(57, 66), (58, 66), (59, 66), (58, 65), (59, 65), (66, 100), (66, 103), (66, 104), (64, 26), (64, 57),
          (18, 103), (19, 104), (18, 99), (42, 48), (39, 50), (44, 49), (4, 20), (1, 28), (2, 28), (3, 20)]
    if round_no >= 5:
        lp += [(39, 49), (12, 39), (17, 30), (5, 20), (56, 67), (59, 67)]
    g.append(('loops', lp))
    return g

def local():
    c = []
    c += helix(25, 45) + helix(76, 82) + helix(84, 106)
    c += [(79, 83), (80, 83), (80, 84), (81, 84), (82, 85), (81, 85)]
    c += [(9, 12), (11, 14), (58, 61), (59, 62), (60, 63), (61, 64), (62, 65), (64, 67), (73, 76), (74, 77), (75, 78), (74, 78),
          (3, 6), (2, 6), (3, 7), (18, 21), (19, 22), (20, 23), (21, 25), (22, 25), (22, 26), (23, 26),
          (45, 48), (46, 49), (47, 50), (47, 51), (48, 51), (104, 107), (105, 108), (103, 107)]
    return c

def all_contacts(round_no=4):
    c = local()
    for _, p in groups(round_no): c += p
    return clean(c)

if __name__ == '__main__':
    rn = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    c = all_contacts(rn)
    print('round', rn, 'n contacts', len(c), 'target', 2 * L, '-', 3 * L, ' with sep>=6:', sum(1 for i, j in c if j - i >= 6))
    X = fold(L, c)
    r = report(X, c)
    segs = dict(Nt=(1, 4), b1=(5, 10), b2=(13, 18), hA=(25, 45), b3=(51, 58), lp=(59, 67), b4=(68, 73), hB=(76, 82), hC=(84, 106))
    seg_dist(r, segs)
    print('close but not predicted:', near_missing(r, c, n=50))
    far = sorted([(i, j, round(float(r[i-1, j-1]), 1)) for i, j in c if r[i-1, j-1] > 7.5 and j - i >= 6], key=lambda t: -t[2])
    print('predicted long-range contacts held only at the limit:', far)
    pickle.dump((c, r), open('r%d.pkl' % rn, 'wb'))
