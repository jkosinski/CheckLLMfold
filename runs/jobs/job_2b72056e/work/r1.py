from lib import *

def base():
    c = []
    # helix 23-34 and short helical turn 56-60: i,i+3 / i,i+4
    for i in range(22, 32): c += [(i, i + 3), (i, i + 4)]
    c += [(32, 35)]
    for i in range(55, 58): c += [(i, i + 3), (i, i + 4)]
    # turns
    c += [(7, 10), (8, 11), (7, 11), (18, 21), (19, 22), (34, 37), (37, 40), (45, 48), (46, 49), (51, 54), (52, 55),
          (60, 63), (62, 65), (72, 75)]
    # hairpin strands 1-7 / 11-17 antiparallel (i + j = 18)
    for i in range(1, 7):
        c += [(i, 18 - i), (i, 17 - i), (i, 19 - i)]
    # strand 1-7 parallel with 64-71 (j = i + 64)
    for i in range(1, 8):
        c += [(i, i + 64), (i, i + 65)]
    # strand 40-45 antiparallel with 66-72 (i + j = 112)
    for i in range(40, 46):
        c += [(i, 112 - i), (i, 113 - i)]
    # short hairpin 43-45 / 48-50
    c += [(44, 49), (43, 50), (44, 50)]
    return c

if __name__ == '__main__':
    c = base()
    # helix packing on strands (hydrophobic faces)
    c += [(3, 26), (5, 26), (3, 23), (15, 26), (17, 23), (13, 30), (15, 30), (5, 30), (27, 43), (30, 43), (30, 41),
          (23, 56), (26, 56), (56, 61), (61, 67), (56, 67), (45, 61), (50, 59), (43, 56)]
    contacts, X, r = run(c, tag='R1')
