"""Score one run. Usage: evaluate.py RUN_DIR REFERENCE(.pdb|.cif) [CHAIN]
Writes RUN_DIR/metrics.json and RUN_DIR/model_final.pdb (CA trace rebuilt from the final contacts with fold())."""
import sys, os, json, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structutil import parse, cb_ca, contact_map, pairs_at_sep, tm_score, seq_of
from foldcode import fold


def last_json(text):
    blocks = re.findall(r'```json\s*(\{.*?\})\s*```', text, re.S)
    if blocks:
        return json.loads(blocks[-1])
    i = text.rfind('{"rounds"')
    return json.loads(text[i:text.rfind('}') + 1])


def precision(pairs, cm, k):
    sel = pairs[:k]
    return (sum(cm[i - 1, j - 1] for i, j in sel) / len(sel)) if sel else float('nan')


def write_ca_pdb(X, seq, path):
    three = {v: k for k, v in {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
             'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
             'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}.items()}
    with open(path, 'w') as f:
        for k, (x, aa) in enumerate(zip(X, seq)):
            f.write('ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n' % (k + 1, three[aa], k + 1, *x))
        f.write('END\n')


def main():
    run, ref = sys.argv[1], sys.argv[2]
    chain = sys.argv[3] if len(sys.argv) > 3 else None
    res = [r for r in parse(ref, chain) if 'CA' in r['atoms']]
    seq = seq_of(res)
    cb, ca = cb_ca(res)
    cm, _ = contact_map(cb)
    L = len(res)
    ans = last_json(open(os.path.join(run, 'final.txt')).read())
    pairs = [(min(p['i'], p['j']), max(p['i'], p['j'])) for p in ans['pairs']]
    pairs = [p for p in pairs if p[1] - p[0] >= 6 and 1 <= p[0] and p[1] <= L]
    seen, ranked = set(), []
    for p in pairs:
        if p not in seen:
            seen.add(p); ranked.append(p)
    final = [tuple(sorted(c)) for c in ans.get('final_contacts', [])]
    final = sorted(set(c for c in final if c[1] - c[0] >= 3 and c[1] <= L))
    nat6 = pairs_at_sep(cm, 6)
    n_pos6 = (L - 6) * (L - 5) // 2
    lr = [p for p in ranked if p[1] - p[0] >= 24]
    m = {'L': L, 'rounds': ans.get('rounds'), 'n_ranked': len(ranked), 'n_final_contacts': len(final),
         'native_contacts_sep3': len(pairs_at_sep(cm, 3)), 'native_contacts_sep6': len(nat6),
         'native_contacts_sep24': len(pairs_at_sep(cm, 24)),
         'chance_precision_sep6': len(nat6) / n_pos6}
    for name, k in (('L5', L // 5), ('L2', L // 2), ('L', L), ('2L', 2 * L)):
        m['P_sep6_top' + name] = precision(ranked, cm, k)
        m['P_sep24_top' + name] = precision(lr, cm, k)
    m['P_sep6_all_ranked'] = precision(ranked, cm, len(ranked))
    if final:
        f6 = [c for c in final if c[1] - c[0] >= 6]
        f24 = [c for c in final if c[1] - c[0] >= 24]
        m['final_P_sep3'] = precision(final, cm, len(final))
        m['final_P_sep6'] = precision(f6, cm, len(f6))
        m['final_P_sep24'] = precision(f24, cm, len(f24)) if f24 else float('nan')
        m['final_recall_sep6'] = sum(cm[i - 1, j - 1] for i, j in f6) / len(nat6)
        for lo, hi, tag in ((3, 5, 'local3_5'), (6, 11, 'short6_11'), (12, 23, 'medium12_23'), (24, 10 ** 6, 'long24')):
            fb = [c for c in final if lo <= c[1] - c[0] <= hi]
            nb = [p for p in pairs_at_sep(cm, lo, hi)]
            m['final_P_' + tag] = precision(fb, cm, len(fb)) if fb else float('nan')
            m['final_R_' + tag] = (sum(cm[i - 1, j - 1] for i, j in fb) / len(nb)) if nb else float('nan')
        X = fold(L, final)
        m['TM_model_vs_ref'] = tm_score(X, ca)
        m['TM_mirror_vs_ref'] = tm_score(X * [1, 1, -1], ca)
        write_ca_pdb(X, seq, os.path.join(run, 'model_final.pdb'))
    Xn = fold(L, pairs_at_sep(cm, 3))
    m['TM_ceiling_native_contacts'] = tm_score(Xn, ca)
    m['TM_ceiling_mirror'] = tm_score(Xn * [1, 1, -1], ca)
    if 'TM_model_vs_ref' in m:
        m['TM_best_hand'] = max(m['TM_model_vs_ref'], m['TM_mirror_vs_ref'])
    m['TM_ceiling_best_hand'] = max(m['TM_ceiling_native_contacts'], m['TM_ceiling_mirror'])
    json.dump(m, open(os.path.join(run, 'metrics.json'), 'w'), indent=1)
    for k, v in m.items():
        print('%-28s %s' % (k, round(v, 3) if isinstance(v, float) else v))


if __name__ == '__main__':
    main()
