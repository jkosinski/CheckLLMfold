import json
import numpy as np

L = 107
cl = json.load(open('c5.json.contacts.json'))
r = np.load('c5.json.r.npy')
T = dict(((a, b), t) for a, b, t in cl)
long = sorted((t, r[a - 1, b - 1], a, b) for (a, b), t in T.items() if b - a >= 6)
ranked = [(a, b) for t, d, a, b in long]
S = set(ranked)
fill = sorted((r[i - 1, j - 1], i, j) for i in range(1, L + 1) for j in range(i + 6, L + 1) if (i, j) not in S)
last = None
for d, i, j in fill:
    if len(ranked) >= 403:
        break
    ranked.append((i, j)); last = d
print(len(ranked), len(set(ranked)), 'predicted long', len(long), 'last fill dist', last)
out = dict(rounds=5, pairs=[dict(i=int(a), j=int(b)) for a, b in ranked],
           final_contacts=[[int(a), int(b)] for a, b, t in cl])
assert all(p['j'] - p['i'] >= 6 for p in out['pairs']) and len(out['final_contacts']) == 317
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(open('answer.json').read())
