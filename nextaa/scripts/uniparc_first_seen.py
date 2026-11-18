"""UniParc first-seen date of each accession's sequence. Usage: uniparc_first_seen.py IN.tsv(accession col) OUT.tsv"""
import sys, csv, urllib.request, urllib.parse, time, os
UA = {'User-Agent': 'python-urllib/3 research-script'}
accs = sorted({r['accession'] for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')})
done = {}
if os.path.exists(sys.argv[2]):
    done = {r['accession']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
with open(sys.argv[2], 'a') as f:
    if not done: f.write('accession\tupi\tfirst_seen\n')
    for a in accs:
        if a in done: continue
        u = 'https://rest.uniprot.org/uniparc/search?' + urllib.parse.urlencode({'query': a, 'format': 'tsv', 'fields': 'upi,first_seen'})
        for t in range(4):
            try:
                lines = urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60).read().decode().strip().split('\n'); break
            except Exception: time.sleep(3 * (t + 1)); lines = ['', '']
        upi, fs = (lines[1].split('\t') + ['', ''])[:2] if len(lines) > 1 else ('', '')
        f.write(f'{a}\t{upi}\t{fs}\n'); f.flush(); time.sleep(0.2)
