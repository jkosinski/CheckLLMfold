"""At each divergence-test site: residue frequencies among the next-closest old relatives (ranks 2..11 by identity),
excluding the closest one. Tells whether the closest relative's residue is a subfamily feature (shared by near
neighbours) or private to that one relative. Usage: analyze_div_neighbors.py DIVDIR NEW.jsonl REL.jsonl"""
import sys, json, csv, collections
d, fnew, frel = sys.argv[1:4]
hits = collections.defaultdict(list)
for l in open(f'{d}/hits.m8'):
    q, t, fid, alen, qs, qe, ts, te, qcov, ev, bits, qaln, taln = l.rstrip('\n').split('\t')
    if float(qcov) >= 0.5: hits[q].append((float(fid), t, int(qs), qaln, taln))
def colmap(qaln, taln, qs):
    m = {}; qi = qs - 1
    for a, b in zip(qaln, taln):
        if a != '-': m[qi] = b; qi += 1
    return m
qs = list(csv.DictReader(open(f'{d}/q_div.tsv'), delimiter='\t'))
relq = [r for r in csv.DictReader(open(f'{d}/q_rel.tsv'), delimiter='\t') if r['condition'] == 'real']
new = {o['qid']: o for o in map(json.loads, open(fnew)) if 'error' not in o and o.get('valid')}
rel = {(o['new_accession'], o['site'], int(o['cut'])): o for o in map(json.loads, open(frel)) if 'error' not in o and o.get('valid') and o['condition'] == 'real'}
out = collections.defaultdict(list)
for a, b in zip(qs, relq):
    n = new.get(a['qid']); r = rel.get((b['new_accession'], b['site'], int(b['cut'])))
    if not (n and r): continue
    k = int(a['cut'])
    nb = sorted([h for h in hits[a['accession']] if h[1] != a['rel']], reverse=True)[:10]
    col = [colmap(h[3], h[4], h[2]).get(k, '-') for h in nb]
    col = [c for c in col if c != '-']
    if not col: continue
    f_rel = col.count(a['rel_res']) / len(col); f_true = col.count(a['true_res']) / len(col)
    out[(a['site'], bool(r['first_correct']))].append((n['pred'][0] == a['rel_res'], n['pred'][0] == a['true_res'], f_rel, f_true, r['pred'][0] == a['rel_res']))
for key in sorted(out):
    v = out[key]; m = len(v)
    print('%-4s Fable right on relative=%-5s n=%3d | on new: says relative %5.1f%%, says new %5.1f%% | 10 next-closest old relatives: relative residue %5.1f%%, new residue %5.1f%%' % (
        key[0], key[1], m, 100 * sum(x[0] for x in v) / m, 100 * sum(x[1] for x in v) / m, 100 * sum(x[2] for x in v) / m, 100 * sum(x[3] for x in v) / m))
