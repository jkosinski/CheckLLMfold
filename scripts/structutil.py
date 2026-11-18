"""Small structure utilities: parse PDB/mmCIF to per-residue CA/CB, contact maps, TM-score (fixed alignment)."""
import numpy as np

THREE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
         'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
         'TYR': 'Y', 'VAL': 'V', 'MSE': 'M'}


def _residues_from_atoms(atoms):
    """atoms: list of (chain, resseq(str incl. icode), resname, atomname, xyz, bfactor). Returns ordered residues."""
    res = {}
    order = []
    for ch, rs, rn, an, xyz, b in atoms:
        key = (ch, rs)
        if key not in res:
            res[key] = {'chain': ch, 'resseq': rs, 'resname': rn, 'atoms': {}, 'b': b}
            order.append(key)
        res[key]['atoms'].setdefault(an, xyz)
    return [res[k] for k in order]


def parse_pdb(path, chain=None, model=1):
    atoms = []
    cur_model = 1
    for line in open(path):
        if line.startswith('MODEL'):
            cur_model = int(line[10:14])
        if line.startswith('ENDMDL') and cur_model == model:
            break
        if cur_model != model or not line.startswith(('ATOM', 'HETATM')):
            continue
        rn = line[17:20].strip()
        if rn not in THREE:
            continue
        alt = line[16]
        if alt not in (' ', 'A'):
            continue
        ch = line[21]
        if chain and ch != chain:
            continue
        atoms.append((ch, line[22:27].strip(), rn, line[12:16].strip(),
                      np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]), float(line[60:66])))
    return _residues_from_atoms(atoms)


def parse_cif(path, chain=None, model=1):
    cols, rows, in_loop = [], [], False
    for line in open(path):
        s = line.strip()
        if s.startswith('_atom_site.'):
            cols.append(s.split('.', 1)[1].split()[0]); in_loop = True; continue
        if in_loop:
            if s.startswith(('ATOM', 'HETATM')):
                rows.append(s.split())
            elif cols and rows and (s.startswith('#') or s.startswith('loop_') or s.startswith('_')):
                break
    ix = {c: k for k, c in enumerate(cols)}
    atoms = []
    for r in rows:
        if 'pdbx_PDB_model_num' in ix and int(r[ix['pdbx_PDB_model_num']]) != model:
            continue
        rn = r[ix['label_comp_id']]
        if rn not in THREE:
            continue
        alt = r[ix['label_alt_id']] if 'label_alt_id' in ix else '.'
        if alt not in ('.', '?', 'A'):
            continue
        ch = r[ix.get('auth_asym_id', ix['label_asym_id'])]
        if chain and ch != chain:
            continue
        rs = r[ix.get('auth_seq_id', ix['label_seq_id'])]
        an = r[ix.get('label_atom_id', ix.get('auth_atom_id'))].strip('"')
        xyz = np.array([float(r[ix['Cartn_x']]), float(r[ix['Cartn_y']]), float(r[ix['Cartn_z']])])
        b = float(r[ix['B_iso_or_equiv']]) if 'B_iso_or_equiv' in ix else 0.0
        atoms.append((ch, rs, rn, an, xyz, b))
    return _residues_from_atoms(atoms)


def parse(path, chain=None):
    return parse_cif(path, chain) if path.endswith('.cif') else parse_pdb(path, chain)


def seq_of(residues):
    return ''.join(THREE[r['resname']] for r in residues)


def cb_ca(residues):
    """Return (CB-or-CA-for-Gly coords, CA coords) arrays [N,3]; residues missing CA are dropped beforehand."""
    cb, ca = [], []
    for r in residues:
        a = r['atoms']
        ca.append(a['CA'])
        if r['resname'] == 'GLY' or 'CB' not in a:
            cb.append(a['CA'])
        else:
            cb.append(a['CB'])
    return np.array(cb), np.array(ca)


def contact_map(cb, cutoff=8.0):
    d = np.linalg.norm(cb[:, None] - cb[None], axis=-1)
    return d < cutoff, d


def pairs_at_sep(cmap, minsep, maxsep=None):
    L = len(cmap)
    out = []
    for i in range(L):
        for j in range(i + minsep, L):
            if maxsep is not None and j - i > maxsep:
                break
            if cmap[i, j]:
                out.append((i + 1, j + 1))
    return out


def kabsch(P, Q):
    """Superpose P onto Q, return rotated+translated P."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return (P - pc) @ R.T + qc


def tm_score(model, ref, Lnorm=None):
    """TM-score of model vs ref (CA coords, same residue order), TM-score program style search over seed fragments."""
    L = len(ref) if Lnorm is None else Lnorm
    d0 = max(1.24 * (L - 15) ** (1 / 3) - 1.8, 0.5) if L > 21 else 0.5
    n = len(ref)
    best = 0.0
    for frac in (1.0, 0.5, 0.25, 0.125):
        Lf = max(int(n * frac), 4)
        step = max(1, Lf // 4) if frac < 1 else n
        for s in range(0, n - Lf + 1, step):
            idx = np.arange(s, s + Lf)
            for _ in range(20):
                if len(idx) < 3:
                    break
                M = kabsch_on(model, ref, idx)
                d = np.linalg.norm(M - ref, axis=1)
                tm = (1 / (1 + (d / d0) ** 2)).sum() / L
                best = max(best, tm)
                cut = d0 + 1.0
                new = np.where(d < cut)[0]
                while len(new) < 3 and cut < 20:
                    cut += 0.5; new = np.where(d < cut)[0]
                if np.array_equal(new, idx):
                    break
                idx = new
    return best


def kabsch_on(P, Q, idx):
    pc, qc = P[idx].mean(0), Q[idx].mean(0)
    H = (P[idx] - pc).T @ (Q[idx] - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return (P - pc) @ R.T + qc
