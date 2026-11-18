"""UniProtKB status of each sequence via UniParc cross-references: 'active' if any UniProtKB entry is active,
'removed' if it was in UniProtKB but all entries are now inactive (date of last removal given), 'never' otherwise.
Usage: uniprot_status.py FASTA ACC1,ACC2,... -> prints TSV; cached in data/uniprot_status_cache.tsv"""
import sys, os, json, urllib.request, urllib.parse, time, hashlib
CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'uniprot_status_cache.tsv')
UA = {'User-Agent': 'python-urllib/3 research-script'}
def status(seq, cache={}):
    if not cache and os.path.exists(CACHE):
        for l in open(CACHE):
            k, *v = l.rstrip('\n').split('\t'); cache[k] = v
    h = hashlib.md5(seq.encode()).hexdigest()
    if h in cache: return cache[h]
    url = 'https://rest.uniprot.org/uniparc/sequence?' + urllib.parse.urlencode({'sequence': seq, 'format': 'json'})
    for t in range(4):
        try:
            d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120)); break
        except urllib.error.HTTPError as e:
            if e.code == 404: d = {}; break
            time.sleep(3 * (t + 1))
        except Exception: time.sleep(3 * (t + 1))
    xr = [x for x in d.get('uniParcCrossReferences', []) if x.get('database', '').startswith('UniProtKB')]
    first = min((x.get('created', '') for x in xr), default='')
    if any(x.get('active') for x in xr): st = 'active'
    elif xr: st = 'removed'
    else: st = 'never'
    last = max((x.get('lastUpdated', '') for x in xr), default='')
    v = [st, first, last, d.get('oldestCrossRefCreated', '')]
    cache[h] = v
    with open(CACHE, 'a') as f: f.write('\t'.join([h] + v) + '\n')
    return v
if __name__ == '__main__':
    seqs = {}; h = None
    for l in open(sys.argv[1]):
        if l.startswith('>'): h = l[1:].strip()
        else: seqs[h] = l.strip()
    print('accession\tuniprot\tuniprot_first\tuniprot_last_update\tuniparc_first')
    for a in sys.argv[2].split(','):
        print(a, *status(seqs[a]), sep='\t', flush=True); time.sleep(0.2)
