from rounds import *
from r1 import base, packing1
from r2 import turns2, packing2

def revise3(cs):
    # dropped: contacts sitting exactly at the 8 A limit that pull unrelated segments together
    for i, j in [(34, 71), (23, 61), (17, 56), (8, 71), (19, 56), (26, 55), (24, 52), (15, 27), (17, 27), (2, 18),
                 (50, 67), (30, 67), (13, 30), (48, 59), (9, 12), (30, 36)]:
        cs.rm(i, j)
    # added: pairs the model puts at 4-6.5 A that fit the side-chain pattern
    for i, j, w in [(21, 55, .45), (19, 55, .35), (21, 54, .35), (1, 19, .4), (17, 25, .4), (17, 24, .35), (15, 25, .3),
                    (39, 72, .4), (39, 74, .3), (50, 57, .35), (46, 59, .3), (55, 61, .3), (24, 55, .3),
                    (3, 15 + 2, .0), (3, 30, .4), (26, 30 + 13, .35), (26, 61, .35), (56, 59 + 2, .0)]:
        if w > 0: cs.add(i, j, w)

if __name__ == '__main__':
    cs = CS(); base(cs); packing1(cs); turns2(cs); packing2(cs); revise3(cs)
    run(cs, 'r3')
