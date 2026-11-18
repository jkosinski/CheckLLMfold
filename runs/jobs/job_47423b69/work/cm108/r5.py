import json
from lib import *
import r4

REMOVE = {(7,42),(3,46),(16,70)}                 # at the limit, weakly motivated, on over-subscribed residues
ADD = {(94,97),(94,98),(22,51),(62,65),(104,108)}  # helix N-cap of 95-107, strand start 22 next to 51/52, helix ends

def contacts():
    s = set(norm(r4.contacts())) - set(norm(REMOVE))
    return s | set(norm(ADD))

# confidence tiers for the ranking
def tiers():
    t = {}
    def put(pairs, v):
        for p in norm(pairs):
            t[p] = max(t.get(p, 0), v)
    put(anti(76, 82, 167, w=0), 10.0)     # hairpin 76-82 / 85-91, tight NG turn
    put(anti(76, 82, 167, w=1), 9.0)
    put(par(23, 28, 30, w=0), 8.5)        # 23-28 / 53-58
    put(par(4, 8, 50, w=0), 8.0)          # 3-8 / 53-58
    put(anti(23, 28, 105, w=0), 8.0)      # 23-28 / 77-82
    put(par(23, 28, 30, w=1), 7.5)
    put(par(3, 8, 50, w=1), 7.0)
    put(anti(23, 28, 105, w=1), 7.0)
    # helix 34-49 on its two flanking strands
    put({(38,27),(41,25),(42,25),(45,23),(46,23),(41,55),(42,55),(45,53),(45,55),(46,53),(49,53),(38,57),(27,41),(28,38),(37,28)}, 6.5)
    # helix 95-107 on the hairpin
    put({(99,91),(99,89),(102,87),(103,87),(103,80),(99,78),(106,85),(86,106),(86,107),(85,107),(95,91),(98,91),(102,80),(106,82)}, 6.0)
    # helix 11-19 on strand 3-8
    put({(12,8),(12,6),(16,4),(16,6),(5,16),(4,17),(5,17),(7,12),(7,13),(6,13)}, 5.5)
    return t

if __name__ == '__main__':
    c, r = run(contacts(), 'r5')
    C = set(c); t = tiers()
    assert all(j - i >= 3 for i, j in c) and 2 * L <= len(c) <= 3 * L
    # score predicted long-range contacts: tier (default 4), small bonus for being comfortably close in the model
    scored = []
    for i, j in c:
        if j - i < 6: continue
        scored.append((t.get((i, j), 4.0) + 0.1 * (8.0 - min(r[i-1, j-1], 8.0)), i, j))
    # fill: unpredicted pairs, by support from neighbouring predicted contacts and model distance
    fill = []
    for i in range(1, L + 1):
        for j in range(i + 6, L + 1):
            if (i, j) in C: continue
            sup = sum(1 for a in range(-2, 3) for b in range(-2, 3)
                      if abs(a) + abs(b) <= 2 and (min(i+a, j+b), max(i+a, j+b)) in C and abs(j + b - i - a) >= 6)
            fill.append((0.25 * sup - 0.25 * r[i-1, j-1], i, j))
    scored.sort(reverse=True); fill.sort(reverse=True)
    N = 406
    ranked = [(i, j) for _, i, j in scored] + [(i, j) for _, i, j in fill]
    ranked = ranked[:N]
    assert len(ranked) == N and len(set(ranked)) == N and all(j - i >= 6 for i, j in ranked)
    print('predicted long-range:', len(scored), '| filled from model:', N - len(scored))
    print('first 25:', ranked[:25])
    print('first fills:', [(i, j, round(float(r[i-1, j-1]), 1)) for _, i, j in fill[:25]])
    print('last fills :', [(i, j, round(float(r[i-1, j-1]), 1)) for _, i, j in fill[N - len(scored) - 10:N - len(scored)]])
    out = {'rounds': 5, 'pairs': [{'i': i, 'j': j} for i, j in ranked], 'final_contacts': [[i, j] for i, j in c]}
    with open('answer.json', 'w') as f: json.dump(out, f, separators=(',', ':'))
    print('final_contacts', len(c), 'pairs', len(ranked))
