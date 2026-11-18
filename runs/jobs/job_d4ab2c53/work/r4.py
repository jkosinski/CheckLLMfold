from rounds import *
from r1 import base, packing1
from r2 import turns2, packing2
from r3 import revise3

def revise4(cs):
    # dropped: weak guesses the model only reaches at the 8 A limit, or that pair an inward with an outward strand side chain
    for i, j in [(24, 55), (13, 26), (59, 62), (39, 74), (46, 59), (15, 25), (27, 52), (34, 41), (55, 61)]:
        cs.rm(i, j)
    # added: loop-to-helix-start and chain-end pairs that the model holds at 4-6.5 A
    for i, j, w in [(22, 56, .45), (22, 54, .35), (1, 22, .35), (41, 73, .4), (38, 73, .3), (47, 59, .3),
                    (3, 25, .3), (27, 56, .35), (53, 59, .3), (36, 72, .3), (15, 29, .3), (1, 21, .3)]:
        cs.add(i, j, w)

def build_all():
    cs = CS(); base(cs); packing1(cs); turns2(cs); packing2(cs); revise3(cs); revise4(cs)
    return cs

if __name__ == '__main__':
    cs = build_all()
    X, r = run(cs, 'r4')
    C = sorted(cs.S)
    assert all(1 <= i < j <= L and j - i >= 3 for i, j in C) and 2 * L <= len(C) <= 3 * L
    lr = [(i, j) for i, j in C if j - i >= 6]
    ranked = sorted(lr, key=lambda p: (-cs.S[p], r[p[0] - 1, p[1] - 1]))
    extra = sorted([(i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1) if (i, j) not in cs.S],
                   key=lambda p: r[p[0] - 1, p[1] - 1])
    N = 286
    ranked = ranked + extra[:N - len(ranked)]
    assert len(ranked) == N and len(set(ranked)) == N
    print('predicted sep>=6: %d, filled from model: %d, max filled distance %.1f' % (len(lr), N - len(lr), r[ranked[-1][0] - 1, ranked[-1][1] - 1]))
    out = {"rounds": 4, "pairs": [{"i": int(i), "j": int(j)} for i, j in ranked], "final_contacts": [[int(i), int(j)] for i, j in C]}
    json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
    print(json.dumps(out, separators=(',', ':')))
