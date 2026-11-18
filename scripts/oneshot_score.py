"""Score the one-shot runs (scripts/oneshot.py) -> results/oneshot/summary.tsv.
contacts: precision of the ranked pairs at top L/5, L/2, L for |i-j| >= 6 and >= 24 (as scripts/evaluate.py).
coords: TM-score of the C-alpha model against the reference (as given, and mirrored), and the same precisions with
pairs ranked by model C-alpha distance. Each call is one row; thinking tokens are kept so calls without thinking can be
compared with the rest."""
import os, sys, json, csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structutil import parse, cb_ca, contact_map, tm_score

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = {'1UBQ': ('data/1UBQ/1ubq.pdb', 'A')}
for line in open(os.path.join(ROOT, 'data/famous.tsv')):
    name, f, ch = line.split()
    REF[name.rsplit('_', 1)[1]] = (f, ch)
for line in open(os.path.join(ROOT, 'data/targets.tsv')):
    a = line.split('\t')[0]
    REF[a] = ('data/%s/%s.pdb' % (a, a), None)


def ref(target):
    f, ch = REF[target]
    res = [r for r in parse(os.path.join(ROOT, f), ch) if 'CA' in r['atoms']]
    cb, ca = cb_ca(res)
    cm, _ = contact_map(cb)
    return cm, ca


def precisions(ranked, cm, L):
    out = {}
    lr = [p for p in ranked if p[1] - p[0] >= 24]
    for name, k in (('L5', L // 5), ('L2', L // 2), ('L', L)):
        for tag, lst in (('sep6', ranked), ('sep24', lr)):
            sel = lst[:k]
            out['P_%s_top%s' % (tag, name)] = sum(cm[i - 1, j - 1] for i, j in sel) / len(sel) if sel else float('nan')
    return out


rows = []
cache = {}
for mode in ('contacts', 'coords'):
    p = os.path.join(ROOT, 'results/oneshot/%s_low.jsonl' % mode)
    if not os.path.exists(p):
        continue
    for l in open(p):
        o = json.loads(l)
        if 'error' in o or o.get('is_error'):
            continue
        t = o['target']
        if t not in cache:
            cache[t] = ref(t)
        cm, ca = cache[t]
        L = len(ca)
        r = dict(target=t, mode=mode, qid=o['qid'], thinking_tokens=o['thinking_tokens'], n_parsed=len(o['parsed']))
        if mode == 'contacts':
            seen, ranked = set(), []
            for i, j in o['parsed']:
                if (i, j) not in seen and j <= L:
                    seen.add((i, j)); ranked.append((i, j))
            r.update(precisions(ranked, cm, L))
        else:
            xyz = {int(n): (x, y, z) for n, x, y, z in o['parsed']}
            if sorted(xyz) != list(range(1, L + 1)):
                r['note'] = 'residues %d of %d' % (len(xyz), L)
                rows.append(r); continue
            X = np.array([xyz[k] for k in range(1, L + 1)])
            r['TM'] = tm_score(X, ca)
            r['TM_mirror'] = tm_score(X * [1, 1, -1], ca)
            d = np.linalg.norm(X[:, None] - X[None], axis=-1)
            ranked = sorted(((i + 1, j + 1) for i in range(L) for j in range(i + 6, L)), key=lambda p: d[p[0] - 1, p[1] - 1])
            r.update(precisions(ranked, cm, L))
            bonds = np.linalg.norm(X[1:] - X[:-1], axis=1)
            r['ca_ca_mean'] = float(bonds.mean())
        rows.append(r)

cols = ['target', 'mode', 'qid', 'thinking_tokens', 'n_parsed', 'TM', 'TM_mirror', 'ca_ca_mean',
        'P_sep6_topL5', 'P_sep6_topL2', 'P_sep6_topL', 'P_sep24_topL5', 'P_sep24_topL2', 'P_sep24_topL', 'note']
with open(os.path.join(ROOT, 'results/oneshot/summary.tsv'), 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', restval='')
    w.writeheader()
    for r in rows:
        w.writerow({k: ('%.3f' % v if isinstance(v, float) else v) for k, v in r.items()})

order = ['1UBQ', '256B', 'A0A916UKZ0', '1TEN', 'A0A8J2YCU1', '2TRX', 'A0A5D0P002']
for mode in ('contacts', 'coords'):
    print('\n' + mode)
    print('%-11s %3s %-20s %-22s %-22s %s' % ('target', 'n', 'thinking tokens', 'top-L |i-j|>=6', 'top-L |i-j|>=24',
                                              'TM (mirror)' if mode == 'coords' else ''))
    for t in order:
        rs = [r for r in rows if r['target'] == t and r['mode'] == mode and 'P_sep6_topL' in r]
        if not rs:
            continue
        f = lambda k: ', '.join('%.2f' % r[k] for r in rs)
        tm = ', '.join('%.2f (%.2f)' % (r['TM'], r['TM_mirror']) for r in rs) if mode == 'coords' else ''
        print('%-11s %3d %-20s %-22s %-22s %s' % (t, len(rs), ','.join(str(r['thinking_tokens']) for r in rs),
                                                  f('P_sep6_topL'), f('P_sep24_topL'), tm))
