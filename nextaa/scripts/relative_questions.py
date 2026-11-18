"""For each divergence-test site, a question on the closest OLD relative itself: prefix = relative up to the aligned
position, target = relative's next 10 residues. Tests whether the relative was memorized.
Usage: relative_questions.py DIVDIR (reads q_div.tsv, hits.m8, old.fasta; writes q_rel.tsv) """
import sys, os, csv, random, zlib
d = sys.argv[1]
old = {}; h = None
for l in open(os.path.join(d, 'old.fasta')):
    if l.startswith('>'): h = l[1:].strip()
    else: old[h] = l.strip()
aln = {}
for l in open(os.path.join(d, 'hits.m8')):
    q, t, fid, alen, qs, qe, ts, te, qcov, ev, bits, qaln, taln = l.rstrip('\n').split('\t')
    aln[(q, t)] = (int(qs), int(ts), qaln, taln)
def shuffled(seq, acc):
    r = list(seq); random.Random(zlib.crc32(acc.encode())).shuffle(r); return ''.join(r)
rows = []
for r in csv.DictReader(open(os.path.join(d, 'q_div.tsv')), delimiter='\t'):
    qs, ts, qaln, taln = aln[(r['accession'], r['rel'])]
    qi, ti, m = qs - 1, ts - 1, {}
    for a, b in zip(qaln, taln):
        if a != '-' and b != '-': m[qi] = ti
        if a != '-': qi += 1
        if b != '-': ti += 1
    k = int(r['cut']); kk = m[k]; rel = old[r['rel']]
    assert rel[kk] == r['rel_res'], (r['qid'], rel[kk], r['rel_res'])
    base = dict(set='relself', accession=r['rel'], pfam=r['pfam'], cut=kk, site=r['site'], new_accession=r['accession'],
                rel_created=r['rel_created'], best_identity=r['best_identity'])
    rows.append(dict(qid=f"relself|{r['rel']}|real|{kk}", condition='real', prefix=rel[:kk], target=rel[kk:kk + 10], **base))
    if r['site'] == 'cons':  # one shuffled control per relative
        s = shuffled(rel, r['rel'])
        rows.append(dict(qid=f"relself|{r['rel']}|shuffled|{kk}", condition='shuffled', prefix=s[:kk], target=s[kk:kk + 10], **base))
keys = list(rows[0].keys())
with open(os.path.join(d, 'q_rel.tsv'), 'w') as f:
    w = csv.DictWriter(f, keys, delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows), 'questions on', len({r['accession'] for r in rows}), 'old relatives')
