import json, sys
import numpy as np
from lib import SEQ, fold, report
from r2 import run, norm

L = len(SEQ)

ADD = {
    3: [  # helix1-helix2 pairs one turn off register that the model placed close, apolar on both sides
        (14, 26), (10, 30), (3, 37), (17, 23), (7, 34),
        # helix1-helix4 neighbours of predicted pairs, model < 5.5 A
        (10, 102), (10, 99), (6, 102), (13, 98), (16, 94), (16, 90),
        # helix3-helix4 near the 82 turn
        (79, 86), (78, 89), (62, 101),
        # loop 43-56 against the helix 2 end / helix 3 start
        (48, 57), (44, 57), (41, 57), (38, 60), (45, 53), (45, 52), (48, 60),
        # helix2-helix3 flanks
        (29, 72), (29, 69), (34, 61)],
    4: [(79, 85), (14, 72), (23, 91), (48, 56), (40, 57), (42, 48), (67, 98), (9, 102)],
    5: [],
}
DROP = {
    3: [  # strained at 8.0 A: one-turn-offset pairs and pairs with a polar partner
        (75, 94), (75, 87), (68, 101), (65, 105), (37, 65), (37, 58), (26, 68), (7, 98), (27, 72), (6, 37), (13, 30), (17, 27),
        # helix2-helix4 pairs the model cannot bring together
        (30, 94), (33, 98)],
    4: [(40, 105), (37, 101), (23, 87), (65, 102), (48, 61), (3, 58)],
    5: [],
}

def make(n):
    c = set(map(tuple, json.load(open('c2.json'))))
    for k in range(3, n + 1):
        c |= set(norm(ADD[k])); c -= set(norm(DROP[k]))
    return sorted(c)

if __name__ == '__main__':
    n = int(sys.argv[1])
    run(n, make(n))
