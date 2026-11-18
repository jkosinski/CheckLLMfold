"""Build the question sets for the recognition probe (scripts/recognize.py).
- data/recognition/contact7.tsv: the seven proteins of the contact test.
- data/recognition/nextaa.tsv: proteins from the next-residue test - the 12 well-known proteins, 40 pre-cutoff proteins
  with no UniRef50 relative, 40 pre-cutoff proteins in UniRef50 clusters of more than 50, and the 56 post-cutoff
  proteins of Tu & Chen's D set."""
import os, csv, random

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NA = os.path.join(ROOT, 'nextaa', 'data')
OUT = os.path.join(ROOT, 'data', 'recognition')
os.makedirs(OUT, exist_ok=True)


def read(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


def write(path, rows, cols):
    with open(path, 'w') as f:
        w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)
    print(path, len(rows))


UBQ = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
c7 = [dict(qid='contact7|1UBQ', set='famous', accession='P0CG48', truth='ubiquitin', organism='Homo sapiens',
           pfam='PF00240', pdb='1UBQ', sequence=UBQ)]
truth = {'cytb562_256B': ('P0ABE7', 'cytochrome b562', 'Escherichia coli', 'PF07361', '256B'),
         'fn3_tenascin_1TEN': ('P24821', 'tenascin FN3 domain 3', 'Homo sapiens', 'PF00041', '1TEN'),
         'thioredoxin_2TRX': ('P0AA25', 'thioredoxin 1', 'Escherichia coli', 'PF00085', '2TRX')}
for r in csv.reader(open(os.path.join(ROOT, 'data', 'famous_targets.tsv')), delimiter='\t'):
    if r[0] in truth:
        a, t, o, p, pdb = truth[r[0]]
        c7.append(dict(qid='contact7|' + pdb, set='famous', accession=a, truth=t, organism=o, pfam=p, pdb=pdb, sequence=r[3]))
duf = {r['Entry']: r for r in read(os.path.join(ROOT, 'data', 'candidates', 'uniprot_duf.tsv'))}
for r in csv.reader(open(os.path.join(ROOT, 'data', 'targets.tsv')), delimiter='\t'):
    d = duf[r[0]]
    c7.append(dict(qid='contact7|' + r[0], set='afdb', accession=r[0], truth=d['Protein names'], organism=d['Organism'],
                   pfam='', pdb='', sequence=r[3]))
cols = ['qid', 'set', 'accession', 'truth', 'organism', 'pfam', 'pdb', 'sequence']
write(os.path.join(OUT, 'contact7.tsv'), c7, cols)

rng = random.Random(1)
u50 = {r['accession']: int(r['member_count']) for r in read(os.path.join(NA, 'main_uniref50.tsv'))}
main = read(os.path.join(NA, 'confirm_sequences.tsv'))
single = sorted([r for r in main if u50.get(r['accession']) == 1], key=lambda r: r['accession'])
big = sorted([r for r in main if u50.get(r['accession'], 0) > 50], key=lambda r: r['accession'])
rows = []
for r in read(os.path.join(NA, 'famous_sequences.tsv')):
    rows.append(dict(qid='famous|' + r['accession'], set='famous', accession=r['accession'], truth=r['name'],
                     organism='', pfam='', uniref50='', sequence=r['sequence']))
for name, pool in (('pre_single', single), ('pre_big', big)):
    for r in rng.sample(pool, 40):
        rows.append(dict(qid='%s|%s' % (name, r['accession']), set=name, accession=r['accession'], truth='',
                         organism=r['organism'], pfam=r['pfam'], uniref50=u50[r['accession']], sequence=r['sequence']))
for r in read(os.path.join(NA, 'd_set_sequences.tsv')):
    rows.append(dict(qid='post_D|' + r['accession'], set='post_D', accession=r['accession'], truth='',
                     organism=r['organism'], pfam=r['pfam'], uniref50='', sequence=r['sequence']))
write(os.path.join(OUT, 'nextaa.tsv'), rows, ['qid', 'set', 'accession', 'truth', 'organism', 'pfam', 'uniref50', 'sequence'])
