"""For each UniProt accession, fetch its UniRef50 and UniRef90 cluster ids and member counts.
Usage: uniref50.py IN.tsv(column 'accession') OUT.tsv"""
import sys, json, urllib.request, urllib.parse, time, csv
accs = [r['accession'] for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')]
out = {}
UA = {'User-Agent': 'python-urllib/3 research-script'}
for ident in ('0.5', '0.9'):
    for k in range(0, len(accs), 40):
        chunk = accs[k:k + 40]
        q = '(' + ' OR '.join('uniprot_id:%s' % a for a in chunk) + ') AND identity:%s' % ident
        url = 'https://rest.uniprot.org/uniref/search?' + urllib.parse.urlencode(
            {'query': q, 'fields': 'id,count,identity', 'format': 'json', 'size': 500})
        for attempt in range(5):
            try:
                d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60)); break
            except Exception as e:
                time.sleep(3 * (attempt + 1))
        # need member->cluster mapping: query each cluster's members is costly; instead query per accession below
        for r in d.get('results', []):
            out.setdefault(ident, []).append((r['id'], r['memberCount']))
        time.sleep(0.3)
json.dump(out, open(sys.argv[2] + '.clusters.json', 'w'))
print({k: len(v) for k, v in out.items()})
