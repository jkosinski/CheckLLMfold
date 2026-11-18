from lib import *
import json
c = [tuple(x) for x in json.load(open('r4.json'))]
add = [(73, 88), (3, 58), (79, 85), (79, 86), (22, 26)]
rem = [(10, 88)]
c = [x for x in c if x not in rem] + add
X, r, c = run(c, 'r5')
json.dump(c, open('r5.json', 'w'))
C = set(c)

H = [(13, 22), (25, 42), (75, 82), (84, 107)]
def same_helix(i, j):
    return any(a <= i <= b and a <= j <= b for a, b in H)
def inseg(i, a, b): return a <= i <= b

def prior(i, j):
    s = i + j
    # hairpin 51-58 / 65-73
    if inseg(i, 50, 59) and inseg(j, 64, 74):
        if s in (123, 124): return 0.85 - 0.004 * (j - i)
        if s in (122, 125): return 0.50 - 0.004 * (j - i)
    # helix 13-22 / helix 25-42 antiparallel pair
    if inseg(i, 11, 22) and inseg(j, 25, 38):
        if s in (47, 48, 49) and (i, j) in C: return 0.55 - 0.003 * (j - i)
        if (i, j) in C: return 0.42 - 0.003 * (j - i)
    # helix 75-82 / start of helix 84-107
    if inseg(i, 74, 82) and inseg(j, 84, 92) and (i, j) in C:
        return 0.50 - 0.012 * (j - i)
    # strand 5-10 / strand 52-57
    if inseg(i, 3, 11) and inseg(j, 51, 59) and (i, j) in C:
        return (0.40 if s == 62 else 0.33)
    # loop 43-50 closing
    if (i, j) in C and inseg(i, 38, 44) and inseg(j, 44, 51): return 0.34
    if (i, j) in C and inseg(i, 70, 72) and inseg(j, 77, 78): return 0.34
    if (i, j) in C: return 0.26
    return None

rows = []
for i in range(1, L + 1):
    for j in range(i + 6, L + 1):
        p = prior(i, j)
        d = r[i - 1, j - 1]
        if p is not None:
            score = p + 0.01 * (8 - min(d, 8)) / 4
        else:
            if same_helix(i, j): continue
            score = 0.20 - 0.02 * d      # filler: unpredicted pairs the model packs closest
        rows.append((score, i, j))
rows.sort(reverse=True)
top = rows[:406]
print('predicted sep>=6:', sum(1 for i, j in c if j - i >= 6), 'ranked:', len(top),
      'of which in final contacts:', sum(1 for s, i, j in top if (i, j) in C))
print('top 40:', [(i, j) for s, i, j in top[:40]])
print('filler cut-off model distance: %.1f' % r[top[-1][1] - 1, top[-1][2] - 1])
out = {"rounds": 5, "pairs": [{"i": i, "j": j} for s, i, j in top], "final_contacts": [[i, j] for i, j in c]}
json.dump(out, open('answer.json', 'w'), separators=(',', ':'))
print(len(c), 'final contacts; all sep>=3:', all(j - i >= 3 for i, j in c), 'unique pairs:', len(set((i, j) for s, i, j in top)))
