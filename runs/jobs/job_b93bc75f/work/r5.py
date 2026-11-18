from lib import *
import pickle, json
L = 89
tier = pickle.load(open('t4.pkl', 'rb'))
# round 5 edits: drop strained / weakly supported pairs, add pairs the model puts close and that fit the packing
for p in [(9, 63), (16, 37), (45, 89), (31, 54)]:
    tier.pop(p, None)
for p in [(70, 83), (67, 86), (20, 37), (32, 54), (24, 33), (26, 32), (67, 87), (44, 66), (33, 53), (74, 82)]:
    tier[p] = 6
c = sorted(tier)
nlong = sum(1 for a, b in c if b - a >= 6)
print(len(c), 'long', nlong)
X = fold(L, c); r = report(X, c)
print(missed(r, c, top=30))
# ranking: predicted long-range contacts by tier (then model distance), then fill with closest unpredicted pairs
pred = sorted([p for p in c if p[1] - p[0] >= 6], key=lambda p: (tier[p], r[p[0] - 1, p[1] - 1]))
S = set(c)
rest = sorted([(i + 1, j + 1) for i in range(L) for j in range(i + 6, L) if (i + 1, j + 1) not in S],
              key=lambda p: r[p[0] - 1, p[1] - 1])
rank = (pred + rest)[:335]
print('ranked', len(rank), 'from predicted', min(len(pred), 335), 'last fill dist %.1f' % r[rank[-1][0] - 1, rank[-1][1] - 1])
assert len(set(rank)) == 335 and all(b - a >= 6 for a, b in rank)
out = {"rounds": 5, "pairs": [{"i": int(a), "j": int(b)} for a, b in rank], "final_contacts": [[int(a), int(b)] for a, b in c]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(json.dumps(out, separators=(',', ':')))
