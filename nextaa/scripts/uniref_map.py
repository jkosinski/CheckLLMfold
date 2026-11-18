"""Map UniProt accessions to their UniRef50 (or UniRef90) cluster and its member count via the UniProt ID-mapping API.
Usage: uniref_map.py IN.tsv(column 'accession') UniRef50|UniRef90 OUT.tsv"""
import sys, json, urllib.request, urllib.parse, time, csv
accs = sorted({r['accession'] for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')})
to = sys.argv[2]
UA = {'User-Agent': 'python-urllib/3 research-script'}
def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers=UA)
    return json.load(urllib.request.urlopen(req, timeout=120))
job = get('https://rest.uniprot.org/idmapping/run',
          urllib.parse.urlencode({'from': 'UniProtKB_AC-ID', 'to': to, 'ids': ','.join(accs)}).encode())['jobId']
while True:
    s = get('https://rest.uniprot.org/idmapping/status/%s' % job)
    if s.get('jobStatus') in (None, 'FINISHED') or 'results' in s: break
    time.sleep(3)
url = 'https://rest.uniprot.org/idmapping/uniref/results/%s?format=json&size=500' % job
rows = []
while url:
    req = urllib.request.Request(url, headers=UA)
    resp = urllib.request.urlopen(req, timeout=120)
    d = json.load(resp)
    for r in d['results']:
        t = r['to']
        rows.append((r['from'], t['id'], t.get('memberCount'), t.get('representativeMember', {}).get('memberId', '')))
    link = resp.headers.get('Link', '')
    url = link.split(';')[0].strip('<> ') if 'rel="next"' in link else None
with open(sys.argv[3], 'w') as f:
    f.write('accession\tcluster\tmember_count\trep\n')
    for r in rows: f.write('\t'.join(map(str, r)) + '\n')
print(len(accs), 'accessions,', len(rows), 'mapped')
