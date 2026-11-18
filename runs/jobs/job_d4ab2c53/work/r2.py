from rounds import *
from r1 import base, packing1

def turns2(cs):
    # local turn / cap contacts in the loops between the strands and helices
    for i, j in [(9, 12), (20, 23), (21, 24), (22, 25), (22, 26), (34, 37), (36, 39), (38, 41), (46, 49),
                 (50, 53), (53, 56), (54, 57), (55, 58), (55, 59), (57, 61), (58, 61), (59, 62), (60, 63), (17, 20)]:
        cs.add(i, j, 0.45)

def packing2(cs):
    for i, j, w in [
        # start of long helix sits over the 13-17 strand (model put 15/16/17 next to 23-26)
        (15, 23, .55), (17, 26, .5), (13, 26, .45), (16, 23, .45), (17, 22, .45), (15, 27, .4), (1, 23, .45), (17, 27, .4),
        # end of long helix over the 41-43 / 69-71 strand ends; I36 in the following loop tucks in there
        (30, 67, .45), (27, 30 + 13, .5), (31, 43, .45), (31, 41, .45), (34, 71, .35), (36, 69, .45), (36, 71, .5),
        (30, 36, .35), (33, 36, .4), (29, 36 + 0, .0),
        # L8 after the first strand next to the 69-71 strand end
        (8, 69, .45), (8, 70, .45), (8, 71, .4),
        # 51-56 loop against start of long helix (K27/E24 vs D52/R54, I23 vs T55/L56)
        (23, 55, .45), (22, 55, .4), (23, 54, .35), (24, 52, .35), (27, 52, .35), (26, 55, .3), (19, 56, .3), (21, 56, .35),
        # short helix 56-60 and I61 on the 43-50 hairpin and 65-67
        (45, 59, .4), (48, 59, .4), (51, 59, .45), (45, 61, .45), (56, 67, .35), (61, 65, .4), (1, 61, .35), (1, 63, .5),
        (43, 56, .4), (50, 67, .35), (45, 50, .4), (40, 73, .5), (39, 73, .35), (41, 72, .5), (2, 18, .3), (1, 18, .4),
    ]:
        if w > 0: cs.add(i, j, w)

if __name__ == '__main__':
    cs = CS(); base(cs); packing1(cs); turns2(cs); packing2(cs)
    run(cs, 'r2')
