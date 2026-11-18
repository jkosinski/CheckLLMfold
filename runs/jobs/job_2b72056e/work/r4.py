import json
from lib import *
from r3 import r3list

def r4list():
    c = set(r3list())
    # V26 cannot reach both I13 and L43; I36 fits L71 better than L69; other pairs held at the 8 A limit
    for p in [(13, 26), (26, 43), (36, 69), (46, 65), (23, 61), (6, 10)]:
        c.discard(p)
    # close in the model and consistent with the strand pairing / turns already predicted
    c |= {(43, 51), (60, 67), (35, 41), (2, 18), (48, 59), (4, 26)}
    return sorted(c)

if __name__ == '__main__':
    contacts, X, r = run(r4list(), tag='R4')
    L = 76
    S = set(contacts)
    # confidence tiers for long-range pairs
    def tier(i, j):
        if 1 <= i <= 7 and 11 <= j <= 17 and i + j == 18: return 0          # hairpin register
        if 1 <= i <= 7 and j == i + 64: return 0                            # parallel register
        if 40 <= i <= 45 and i + j == 112: return 0                         # antiparallel register
        if 1 <= i <= 7 and 11 <= j <= 17 and abs(i + j - 18) == 1: return 1
        if 1 <= i <= 7 and abs(j - i - 64) == 1: return 1
        if 40 <= i <= 46 and 65 <= j <= 73 and abs(i + j - 112) == 1: return 1
        if 42 <= i <= 45 and 48 <= j <= 51: return 2
        hyd = set([1, 3, 5, 13, 15, 17, 23, 26, 27, 30, 36, 41, 43, 45, 50, 56, 59, 61, 67, 69, 71])
        if i in hyd and j in hyd: return 3
        return 4
    longp = [(i, j) for i, j in contacts if j - i >= 6]
    longp.sort(key=lambda p: (tier(*p), r[p[0] - 1, p[1] - 1]))
    rest = sorted([(r[i - 1, j - 1], i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1) if (i, j) not in S])
    pairs = longp + [(i, j) for _, i, j in rest]
    pairs = pairs[:286]
    assert len(set(pairs)) == 286 and all(j - i >= 6 for i, j in pairs)
    print('long predicted', len(longp), 'filled from model', 286 - len(longp))
    out = {"rounds": 4, "pairs": [{"i": int(i), "j": int(j)} for i, j in pairs],
           "final_contacts": [[int(i), int(j)] for i, j in contacts]}
    json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
    print(json.dumps(out, separators=(',', ':')))
