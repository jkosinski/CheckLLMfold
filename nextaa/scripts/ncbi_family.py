"""Download every RefSeq protein whose title contains a family name, with its NCBI creation date.
Usage: ncbi_family.py NAME OUTDIR  -> OUTDIR/NAME.fasta, OUTDIR/NAME.dates.tsv (accession, createdate, length, title)"""
import sys, os, json, time, urllib.request, urllib.parse
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
UA = {'User-Agent': 'python-urllib/3 research-script'}
def get(ep, **p):
    p['tool'] = 'research-script'
    for a in range(6):
        try:
            time.sleep(0.4)
            return urllib.request.urlopen(urllib.request.Request(E + ep, data=urllib.parse.urlencode(p).encode(), headers=UA), timeout=180).read().decode()
        except Exception as e:
            time.sleep(5 * (a + 1)); err = e
    raise err
name, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
fa, dt = os.path.join(out, name + '.fasta'), os.path.join(out, name + '.dates.tsv')
if os.path.exists(dt) and os.path.getsize(dt) > 0: print(name, 'cached'); sys.exit()
r = json.loads(get('esearch.fcgi', db='protein', term=f'{name}[Title] AND refseq[filter]', usehistory='y', retmode='json', retmax=0))['esearchresult']
n, wenv, qk = int(r['count']), r['webenv'], r['querykey']
seqs, rows = [], []
B = 500
for s in range(0, n, B):
    j = json.loads(get('esummary.fcgi', db='protein', WebEnv=wenv, query_key=qk, retstart=s, retmax=B, retmode='json'))['result']
    for uid in j.get('uids', []):
        d = j[uid]; rows.append((d['caption'] + ('.' + str(d.get('accessionversion', '')).split('.')[-1] if False else ''), d['createdate'], d['slen'], d['title']))
    seqs.append(get('efetch.fcgi', db='protein', WebEnv=wenv, query_key=qk, retstart=s, retmax=B, rettype='fasta', retmode='text'))
open(fa, 'w').write(''.join(seqs))
with open(dt, 'w') as f:
    f.write('accession\tcreatedate\tlength\ttitle\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
print(name, n, 'records', len(rows), 'dates', sum(x.count('>') for x in seqs), 'seqs', flush=True)
