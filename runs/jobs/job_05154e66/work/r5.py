from common import *
from r4 import contacts_r4

c = clean(contacts_r4())
# dropped: cross contacts to strand 92-96 that the round-4 model could not satisfy
drop = {(12, 94), (79, 94), (14, 92), (16, 79), (31, 103), (78, 107)}
c = [p for p in c if p not in drop]
# added: pairs the round-4 model placed close and that fit the side-chain facing of the strands / helix start
c += [(7, 72), (30, 72), (8, 35), (13, 62), (15, 60), (11, 62), (9, 66), (99, 106), (98, 107), (75, 105), (65, 72),
      (7, 36), (24, 95), (22, 58)]
X, r, c = run(c, 'r5')

# ---------- ranking of pairs with |i - j| >= 6 ----------
score = {}


def put(i, j, s):
    i, j = min(i, j), max(i, j)
    if j - i >= 6 and 1 <= i and j <= L:
        score[(i, j)] = max(score.get((i, j), -1e9), s)


def strandpair(s, lo, hi, base):
    mid = (lo + hi) / 2.0
    for i in range(lo, hi + 1):
        edge = abs(i - mid) / max(1.0, (hi - lo) / 2.0)   # 0 centre .. 1 end
        for dj, pen in ((0, 0), (1, 8), (-1, 8), (2, 28), (-2, 28)):
            put(i, s - i + dj, base - pen - 4 * edge)


strandpair(37, 6, 16, 100)     # hairpin 6-16 / 21-31
strandpair(175, 76, 85, 95)    # hairpin 76-85 / 90-99
strandpair(73, 7, 16, 90)      # 7-16 / 57-66
strandpair(142, 59, 67, 88)    # 59-67 / 75-83
cs = set(c)
for (i, j) in c:
    if j - i >= 6 and (i, j) not in score:
        if 32 <= i <= 52 or 32 <= j <= 52:
            put(i, j, 62)          # helix against strands
        elif 23 <= i <= 29 and 91 <= j <= 99:
            put(i, j, 50 if j - i == 69 else 44)
        else:
            put(i, j, 47)
# fill with pairs the final model puts close
for i in range(1, L + 1):
    for j in range(i + 6, L + 1):
        if (i, j) not in score:
            put(i, j, 30 - r[i - 1, j - 1])
rank = sorted(score, key=lambda p: -score[p])[:403]
print('n ranked', len(rank), 'of which in final contacts', sum(1 for p in rank if p in cs))
print('top 25', rank[:25])
print('last 10', rank[-10:], [round(score[p], 1) for p in rank[-10:]])
out = {"rounds": 5, "pairs": [{"i": int(i), "j": int(j)} for i, j in rank], "final_contacts": [[int(i), int(j)] for i, j in c]}
assert len(set(rank)) == 403 and all(j - i >= 6 for i, j in rank)
assert 2 * L <= len(c) <= 3 * L, len(c)
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print('final contacts', len(c))
