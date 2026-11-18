"""For each Foldseek result JSON, compute TM-score (normalised by query length) of every hit using Foldseek's own
residue alignment and the returned target CA trace. Usage: foldseek_tm.py MODELDIR RESULT.json [...]"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from structutil import parse, cb_ca, tm_score


def hits_of(d):
    out = []
    for db in d.get('results', []):
        for aln in db.get('alignments') or []:
            for a in (aln if isinstance(aln, list) else [aln]):
                out.append(a)
    return out


def aln_tm(a, qca):
    tca = np.array([float(x) for x in a['tCa'].split(',')]).reshape(-1, 3)
    qi, ti = a['qStartPos'] - 1, a['dbStartPos'] - 1
    qs, ts = [], []
    for cq, ct in zip(a['qAln'], a['dbAln']):
        if cq != '-' and ct != '-':
            qs.append(qi); ts.append(ti)
        if cq != '-':
            qi += 1
        if ct != '-':
            ti += 1
    if len(qs) < 5:
        return 0.0
    return tm_score(tca[ts], qca[qs], Lnorm=len(qca))


def main():
    mdir = sys.argv[1]
    for f in sys.argv[2:]:
        acc = os.path.basename(f).split('.')[0]
        _, qca = cb_ca(parse(os.path.join(mdir, acc + '.pdb')))
        hs = hits_of(json.load(open(f)))
        best = (0, None)
        for a in hs:
            tm = aln_tm(a, qca)
            if tm > best[0]:
                best = (tm, a)
        sig = [a for a in hs if float(a['eval']) < 0.1]
        b = best[1] or {}
        print('%s nhits=%d n(E<0.1)=%d minE=%.3g maxprob=%.2f bestTM=%.2f (%s, E=%s)' % (
            acc, len(hs), len(sig), min([float(a['eval']) for a in hs] or [99]), max([a['prob'] for a in hs] or [0]),
            best[0], b.get('target', '-')[:50], b.get('eval')))


if __name__ == '__main__':
    main()
