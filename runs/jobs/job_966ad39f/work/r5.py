from lib import *
from r4 import groups, local
import pickle, json

DROP = set([(71, 81), (15, 36), (36, 53)])   # polar / strained at the 8 A limit in rounds 3-4

def final_groups():
    return [(n, [q for q in clean(p) if q not in DROP]) for n, p in groups(5)]

def final_contacts():
    c = local()
    for _, p in final_groups(): c += p
    return clean(c)

if __name__ == '__main__':
    c = final_contacts()
    print('round 5 n contacts', len(c), 'target', 2 * L, '-', 3 * L, ' with sep>=6:', sum(1 for i, j in c if j - i >= 6))
    X = fold(L, c)
    r = report(X, c)
    segs = dict(Nt=(1, 4), b1=(5, 10), b2=(13, 18), hA=(25, 45), b3=(51, 58), lp=(59, 67), b4=(68, 73), hB=(76, 82), hC=(84, 106))
    seg_dist(r, segs)
    print('close but not predicted:', near_missing(r, c, n=40))
    pickle.dump((c, r), open('r5.pkl', 'wb'))

    # ranking: confidence groups first (in order), then remaining pairs by model distance,
    # favouring pairs that neighbour a predicted long-range contact (register hedges)
    ranked = []; seen = set()
    for name, p in final_groups():
        for q in p:
            if q[1] - q[0] >= 6 and q not in seen:
                ranked.append(q); seen.add(q)
    print('predicted pairs with sep>=6:', len(ranked))
    LR = [q for q in c if q[1] - q[0] >= 6]
    cand = []
    for i in range(1, L + 1):
        for j in range(i + 6, L + 1):
            if (i, j) in seen: continue
            nb = min(abs(i - a) + abs(j - b) for a, b in LR)
            score = r[i - 1, j - 1] + 0.8 * min(nb, 4)
            cand.append((score, i, j))
    cand.sort()
    for s, i, j in cand:
        if len(ranked) >= 406: break
        ranked.append((i, j))
    print('ranked', len(ranked), 'last fill score', round(float(s), 2))
    assert len(set(ranked)) == 406 and all(j - i >= 6 for i, j in ranked)
    assert all(j - i >= 3 for i, j in c)
    out = dict(rounds=5, pairs=[dict(i=int(i), j=int(j)) for i, j in ranked], final_contacts=[[int(i), int(j)] for i, j in c])
    json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
