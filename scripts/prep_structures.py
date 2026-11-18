"""For each target, take the run shown in the X-thread contact maps (best long-range precision; the only run for single-run
targets), rebuild its CA model from the final contacts, pick the better mirror image, superpose it on the reference with the
TM-score superposition, and write CA-only PDBs: results/x_thread/structures/<target>_ref.pdb and <target>_model.pdb.
Also writes structures.json with TM-scores and the chosen job."""
import os, sys, csv, json
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
from structutil import parse, cb_ca, seq_of
from evaluate import write_ca_pdb, last_json
from foldcode import fold

OUT = os.path.join(ROOT, 'results/x_thread/structures')
os.makedirs(OUT, exist_ok=True)
REFS = {'ubiquitin_1UBQ': ('data/1UBQ/1ubq.pdb', 'A'), 'cytb562_256B': ('data/256B/256b.pdb', 'A'),
        'fn3_tenascin_1TEN': ('data/1TEN/1ten.pdb', 'A'), 'thioredoxin_2TRX': ('data/2TRX/2trx.pdb', 'A')}
for a in ('A0A916UKZ0', 'A0A8J2YCU1', 'A0A5D0P002'):
    REFS[a] = ('data/%s/%s.pdb' % (a, a), None)
rows = [r for r in csv.DictReader(open(os.path.join(ROOT, 'results/summary.tsv')), delimiter='\t') if r['job'].startswith('job_')]


def tm_superpose(model, ref):
    """Return (TM, model superposed on ref) using the TM-score search over seed fragments."""
    L = len(ref)
    d0 = max(1.24 * (L - 15) ** (1 / 3) - 1.8, 0.5)
    best, bestX = -1, None
    for frac in (1.0, 0.5, 0.25, 0.125):
        Lf = max(int(L * frac), 4)
        step = max(1, Lf // 4) if frac < 1 else L
        for s in range(0, L - Lf + 1, step):
            idx = np.arange(s, s + Lf)
            for _ in range(20):
                pc, qc = model[idx].mean(0), ref[idx].mean(0)
                U, S, Vt = np.linalg.svd((model[idx] - pc).T @ (ref[idx] - qc))
                R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
                X = (model - pc) @ R.T + qc
                d = np.linalg.norm(X - ref, axis=1)
                tm = (1 / (1 + (d / d0) ** 2)).sum() / L
                if tm > best:
                    best, bestX = tm, X
                new = np.where(d < d0 + 1)[0]
                if len(new) < 3 or np.array_equal(new, idx):
                    break
                idx = new
    return best, bestX


info = {}
for t, (f, ch) in REFS.items():
    cand = [r for r in rows if r['target'] == t]
    best = max(cand, key=lambda r: float(r['P_sep24_topL']))
    res = [r for r in parse(os.path.join(ROOT, f), ch) if 'CA' in r['atoms']]
    _, ca = cb_ca(res)
    seq = seq_of(res)
    ans = last_json(open(os.path.join(ROOT, 'runs/jobs', best['job'], 'final.txt')).read())
    final = sorted(set(tuple(sorted(c)) for c in ans['final_contacts'] if abs(c[1] - c[0]) >= 3 and max(c) <= len(res)))
    X = fold(len(res), final)
    tm1, X1 = tm_superpose(X, ca)
    tm2, X2 = tm_superpose(X * [1, 1, -1], ca)
    tm, Xs, hand = (tm1, X1, 'as built') if tm1 >= tm2 else (tm2, X2, 'mirrored')
    ca_c = ca - ca.mean(0)
    write_ca_pdb(ca_c, seq, os.path.join(OUT, '%s_ref.pdb' % t))
    write_ca_pdb(Xs - ca.mean(0), seq, os.path.join(OUT, '%s_model.pdb' % t))
    info[t] = {'job': best['job'], 'rep': best['rep'], 'n_runs': len(cand), 'TM': round(tm, 3), 'hand': hand, 'L': len(res)}
    print(t, info[t])
json.dump(info, open(os.path.join(OUT, 'structures.json'), 'w'), indent=1)
