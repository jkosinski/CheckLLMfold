"""Divergence-site test set.
New proteins = RefSeq family members created on/after NEW_FROM; old = created before OLD_BEFORE (likely in training).
For each new protein, MMseqs2 finds old relatives; keep those whose closest old relative is LO-HI identical over >=80%
of the new protein. Cut sites: 'div' where the new residue differs from the closest relative's aligned residue,
'cons' where they agree; the alignment must be gap-free 5 residues either side.
Per site, the residue profile of all other old relatives at that column (excluding the closest one) is recorded.

Usage: divergence_set.py FAMDIR OUTDIR [--lo 0.70 --hi 0.95 --per-family 5 --div 3 --cons 1]"""
import os, sys, csv, glob, random, subprocess, argparse, collections, re, zlib
STD = set('ACDEFGHIKLMNPQRSTVWY')
MMSEQS = os.path.expanduser('~/miniforge3/envs/mmseqs/bin/mmseqs')

def read_fasta(p):
    out = {}; h = None
    for l in open(p):
        if l.startswith('>'): h = l[1:].split()[0].split('.')[0]; out[h] = []
        elif h: out[h].append(l.strip())
    return {k: ''.join(v) for k, v in out.items()}

UNIPARC_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'uniparc_seq_cache.tsv')
_cache = None
def uniparc_oldest(seq):
    """Oldest cross-reference date of this exact sequence in UniParc ('' if UniParc does not have it)."""
    import urllib.request, urllib.parse, json, hashlib, time
    global _cache
    if _cache is None:
        _cache = {}
        if os.path.exists(UNIPARC_CACHE):
            for l in open(UNIPARC_CACHE):
                k, v = l.rstrip('\n').split('\t'); _cache[k] = v
    h = hashlib.md5(seq.encode()).hexdigest()
    if h in _cache: return _cache[h]
    url = 'https://rest.uniprot.org/uniparc/sequence?' + urllib.parse.urlencode({'sequence': seq, 'format': 'json', 'fields': 'upi,first_seen'})
    v = None
    for t in range(4):
        try:
            v = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'python-urllib/3 research-script'}), timeout=60)).get('oldestCrossRefCreated', ''); break
        except urllib.error.HTTPError as e:
            if e.code in (400, 404): v = ''; break
            time.sleep(3 * (t + 1))
        except Exception:
            time.sleep(3 * (t + 1))
    if v is None: v = 'ERROR'
    _cache[h] = v
    with open(UNIPARC_CACHE, 'a') as f: f.write(f'{h}\t{v}\n')
    return v

def colmap(qaln, taln, qstart):
    """query index (0-based) -> (target residue or '-', gapfree flag computed later)"""
    m = {}; qi = qstart - 1
    for a, b in zip(qaln, taln):
        if a != '-':
            m[qi] = b; qi += 1
    return m

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('famdir'); ap.add_argument('out')
    ap.add_argument('--new-from', default='2026/08/01'); ap.add_argument('--old-before', default='2026/06/01')
    ap.add_argument('--lo', type=float, default=0.70); ap.add_argument('--hi', type=float, default=0.95)
    ap.add_argument('--per-family', type=int, default=5); ap.add_argument('--div', type=int, default=3); ap.add_argument('--cons', type=int, default=1)
    ap.add_argument('--old-dir', default='')  # take OLD relatives from this dir (e.g. UniProtKB members) instead of famdir
    ap.add_argument('--extra-queries', default='')  # tsv with accession, sequence (e.g. the D set) searched as well
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    old, new, fam_of, created = {}, {}, {}, {}
    for dt in glob.glob(os.path.join(a.famdir, '*.dates.tsv')):
        fam = os.path.basename(dt).split('.')[0]
        seqs = read_fasta(dt.replace('.dates.tsv', '.fasta'))
        for r in csv.DictReader(open(dt), delimiter='\t'):
            acc = r['accession'].split('.')[0]; s = seqs.get(acc)
            if not s or not (120 <= len(s) <= 800) or set(s) - STD: continue
            created[acc] = r['createdate']; fam_of[acc] = fam
            if r['createdate'] < a.old_before and not a.old_dir: old[acc] = s
            elif r['createdate'] >= a.new_from: new[acc] = s
    if a.old_dir:
        for dt in glob.glob(os.path.join(a.old_dir, '*.dates.tsv')):
            fam = os.path.basename(dt).split('.')[0]
            if not os.path.exists(os.path.join(a.famdir, fam + '.dates.tsv')): continue
            seqs = read_fasta(dt.replace('.dates.tsv', '.fasta'))
            for r in csv.DictReader(open(dt), delimiter='\t'):
                acc = r['accession'].split('.')[0]; s = seqs.get(acc)
                if not s or set(s) - STD: continue
                if r['createdate'] < a.old_before: old[acc] = s; created[acc] = r['createdate']; fam_of[acc] = fam
    extra = {}
    if a.extra_queries:
        for r in csv.DictReader(open(a.extra_queries), delimiter='\t'): extra[r['accession']] = r['sequence']
    print('old', len(old), 'new', len(new), 'extra', len(extra), flush=True)
    with open(os.path.join(a.out, 'old.fasta'), 'w') as f:
        for k, s in old.items(): f.write(f'>{k}\n{s}\n')
    with open(os.path.join(a.out, 'new.fasta'), 'w') as f:
        for k, s in {**new, **extra}.items(): f.write(f'>{k}\n{s}\n')
    res = os.path.join(a.out, 'hits.m8')
    if not os.path.exists(res):
        subprocess.run([MMSEQS, 'easy-search', os.path.join(a.out, 'new.fasta'), os.path.join(a.out, 'old.fasta'), res,
                        os.path.join(a.out, 'tmp'), '-s', '7.5', '--max-seqs', '1000', '-e', '1e-3', '--threads', '6',
                        '--format-output', 'query,target,fident,alnlen,qstart,qend,tstart,tend,qcov,evalue,bits,qaln,taln',
                        '-v', '1'], check=True)
    hits = collections.defaultdict(list)
    for l in open(res):
        q, t, fid, alen, qs, qe, ts, te, qcov, ev, bits, qaln, taln = l.rstrip('\n').split('\t')
        hits[q].append(dict(t=t, fid=float(fid), qs=int(qs), qe=int(qe), qcov=float(qcov), qaln=qaln, taln=taln))
    # best old relative per query (qcov >= 0.8, highest identity)
    summ = []
    for q, hs in hits.items():
        good = [h for h in hs if h['qcov'] >= 0.8]
        best = max(good, key=lambda h: h['fid']) if good else None
        best_any = max(hs, key=lambda h: h['fid'])
        summ.append(dict(accession=q, family=fam_of.get(q, 'extra'), created=created.get(q, ''), length=len(new.get(q) or extra.get(q)),
                         n_old_hits=len(hs), best_id_cov80=best['fid'] if best else 0, best_rel=best['t'] if best else '',
                         best_rel_created=created.get(best['t'], '') if best else '', best_id_any=best_any['fid']))
    for q in list(new) + list(extra):
        if q not in hits: summ.append(dict(accession=q, family=fam_of.get(q, 'extra'), created=created.get(q, ''), length=len(new.get(q) or extra.get(q)), n_old_hits=0, best_id_cov80=0, best_rel='', best_rel_created='', best_id_any=0))
    with open(os.path.join(a.out, 'nearest_old.tsv'), 'w') as f:
        w = csv.DictWriter(f, list(summ[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(summ)
    # site selection for new proteins with a close relative
    cand = [s for s in summ if s['accession'] in new and a.lo <= s['best_id_cov80'] <= a.hi]
    byfam = collections.defaultdict(list)
    for s in cand: byfam[s['family']].append(s)
    rows = []
    for fam, ss in sorted(byfam.items()):
        rng = random.Random(zlib.crc32(fam.encode())); ss.sort(key=lambda s: s['accession']); rng.shuffle(ss); taken = 0; seen_rel = set()
        for s in ss:
            if taken >= a.per_family: break
            if s['best_rel'] in seen_rel: continue  # one new protein per old relative
            q = s['accession']; seq = new[q]
            first = uniparc_oldest(seq)
            if first == 'ERROR' or (first and first < '2026-07-01'): continue  # sequence already public before the cutoff
            s['uniparc_first'] = first or 'not in UniParc'
            hs = [h for h in hits[q] if h['qcov'] >= 0.5]
            best = next(h for h in hs if h['t'] == s['best_rel'])
            bm = colmap(best['qaln'], best['taln'], best['qs'])
            others = [colmap(h['qaln'], h['taln'], h['qs']) for h in hs if h['t'] != best['t']]
            def gapfree(k):
                return all(bm.get(i, '-') != '-' for i in range(k - 5, k + 6)) and all(i in bm for i in range(k - 5, k + 6))
            ks = [k for k in range(30, len(seq) - 10) if k in bm and gapfree(k)]
            div = [k for k in ks if bm[k] != seq[k]]; cons = [k for k in ks if bm[k] == seq[k]]
            if len(div) < a.div or len(cons) < a.cons: continue
            pick = [(k, 'div') for k in sorted(rng.sample(div, a.div))] + [(k, 'cons') for k in sorted(rng.sample(cons, a.cons))]
            for k, typ in pick:
                col = collections.Counter(m[k] for m in others if k in m and m[k] != '-')
                tot = sum(col.values())
                rows.append(dict(qid=f'div|{q}|real|{k}', set='div', accession=q, pfam=fam, condition='real', cut=k,
                                 prefix=seq[:k], target=seq[k:k + 10], site=typ, true_res=seq[k], rel_res=bm[k],
                                 rel=best['t'], rel_created=created.get(best['t'], ''), best_identity=round(s['best_id_cov80'], 3),
                                 n_profile=tot, f_true=round(col[seq[k]] / tot, 3) if tot else '',
                                 f_rel=round(col[bm[k]] / tot, 3) if tot else '',
                                 majority=col.most_common(1)[0][0] if tot else '', created=s['created'], uniparc_first=s['uniparc_first']))
            taken += 1; seen_rel.add(s['best_rel'])
    with open(os.path.join(a.out, 'q_div.tsv'), 'w') as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
    print('candidates', len(cand), 'proteins used', len({r['accession'] for r in rows}), 'questions', len(rows))

if __name__ == '__main__':
    main()
