from lib import *
import pickle, json
T = pickle.load(open('T4.pkl', 'rb'))
def add(pairs, w):
    for i, j in pairs:
        if i > j: i, j = j, i
        if j - i >= 3 and i >= 1 and j <= L: T[(i, j)] = max(T.get((i, j), 0), w)
# pairs the round-4 model put close and that fit the secondary structure
add([(17,90),(18,89),(6,32),(21,51),(22,51),(23,47),(27,43),(35,64),(42,60),(16,86)], 0.45)
add([(66,73),(67,73),(68,74)], 0.45)
# drop the weakest tail guesses that stayed strained
for p in [(91,107),(92,105),(93,104),(4,37)]:
    T.pop(p, None)
c = sorted(T); print(len(c))
X = fold(L, c); r = report(X, c)
print(missed(r, c, n=30))
segs = dict(b1=(6,17), b2=(22,31), h1=(32,53), b3=(56,66), lp=(67,75), h2=(76,87), b4=(91,97), ct=(98,107), nt=(1,5))
segdist(r, segs)

# ranking: my confidence first, model distance as tie-break; then fill from the model
lr = [(p, w) for p, w in T.items() if p[1] - p[0] >= 6]
lr.sort(key=lambda t: (-t[1], r[t[0][0] - 1, t[0][1] - 1]))
ranked = [p for p, _ in lr]
print('own long-range pairs', len(ranked))
def in_helix(i, j):
    return any(a <= i and j <= b for a, b in [(32, 53), (76, 87)])
rest = [((i + 1, j + 1), r[i, j]) for i in range(L) for j in range(i + 6, L)
        if (i + 1, j + 1) not in T and not in_helix(i + 1, j + 1)]
rest.sort(key=lambda t: t[1])
for p, d in rest:
    if len(ranked) >= 403: break
    ranked.append(p)
print('fill cut-off distance', d, 'total', len(ranked), 'unique', len(set(ranked)))
out = {"rounds": 5, "pairs": [{"i": int(i), "j": int(j)} for i, j in ranked],
       "final_contacts": [[int(i), int(j)] for i, j in c]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(len(out['pairs']), len(out['final_contacts']), min(j - i for i, j in ranked))
