"""Build question files (qid, set, accession, pfam, condition, cut, prefix, target, + metadata) for probe.py.
Shuffle and cut conventions are Tu & Chen's: prefix = seq[:cut], target = seq[cut:cut+10];
shuffled = residues shuffled by random.Random(zlib.crc32(accession)), same cut."""
import csv, random, zlib, gzip, json, sys

def shuffled(seq, acc):
    r = list(seq); random.Random(zlib.crc32(acc.encode())).shuffle(r); return ''.join(r)

def q(set_, acc, pfam, cond, seq, cut, **meta):
    s = seq if cond == 'real' else shuffled(seq, acc)
    return dict(qid=f'{set_}|{acc}|{cond}|{cut}', set=set_, accession=acc, pfam=pfam, condition=cond, cut=cut,
                prefix=s[:cut], target=s[cut:cut + 10], **meta)

def write(rows, path):
    keys = list(rows[0].keys())
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(path, 'w') as f:
        w = csv.DictWriter(f, keys, delimiter='\t', restval=''); w.writeheader(); w.writerows(rows)
    print(path, len(rows))

if __name__ == '__main__':
    # check our shuffle reproduces theirs
    main = {r['accession']: r for r in csv.DictReader(open('data/confirm_sequences.tsv'), delimiter='\t')}
    n = ok = 0
    for l in gzip.open('data/tu_outputs/confirm/claude-fable-5-2.jsonl.gz', 'rt'):
        o = json.loads(l); r = main[o['accession']]
        s = r['sequence'] if o['condition'] == 'real' else shuffled(r['sequence'], o['accession'])
        n += 1; ok += s[o['cut']:o['cut'] + 10] == o['target']
        if n >= 400: break
    print('targets reproduced', ok, 'of', n)
    ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open('data/main_uniref50.tsv'), delimiter='\t')}
    rows = []
    # A. well-known proteins: 5 cuts in [30, L-10), seeded by accession; real + shuffled
    for r in csv.DictReader(open('data/famous_sequences.tsv'), delimiter='\t'):
        L = int(r['length']); hi = 76 - 10 if r['accession'] == 'P0CG48' else L - 10
        cuts = sorted(random.Random(zlib.crc32(('W|' + r['accession']).encode())).sample(range(30, hi), 5))
        for c in cuts:
            for cond in ('real', 'shuffled'):
                rows.append(q('famous', r['accession'], '', cond, r['sequence'], c, name=r['name']))
    write(rows, 'data/q_famous.tsv')
    # B. main-set proteins with no UniRef50 relative (pre-cutoff, 116 proteins): their first 2 cuts real, first cut shuffled
    rows = []
    for acc, r in main.items():
        if ur50[acc] != 1: continue
        cuts = list(map(int, r['cuts'].split(',')))
        for c in cuts[:2]: rows.append(q('pre_single', acc, r['pfam'], 'real', r['sequence'], c))
        rows.append(q('pre_single', acc, r['pfam'], 'shuffled', r['sequence'], cuts[0]))
    write(rows, 'data/q_pre_single.tsv')
    # C. D set (post-cutoff, 56 proteins): all 5 cuts real, first 2 cuts shuffled
    rows = []
    for r in csv.DictReader(open('data/d_set_sequences.tsv'), delimiter='\t'):
        cuts = list(map(int, r['cuts'].split(',')))
        for c in cuts: rows.append(q('post_D', r['accession'], r['pfam'], 'real', r['sequence'], c, best_identity=r['best_identity']))
        for c in cuts[:2]: rows.append(q('post_D', r['accession'], r['pfam'], 'shuffled', r['sequence'], c, best_identity=r['best_identity']))
    write(rows, 'data/q_post_D.tsv')
