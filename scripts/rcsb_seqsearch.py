"""RCSB MMseqs2 sequence search against the PDB (experimental + computed off).
Usage: rcsb_seqsearch.py IN.tsv(acc<TAB>seq) OUT.tsv ; writes acc, n_hits, best hit id, evalue, identity."""
import sys, json, urllib.request, time

URL = 'https://search.rcsb.org/rcsbsearch/v2/query'


def search(seq, evalue=1.0):
    q = {"query": {"type": "terminal", "service": "sequence",
                   "parameters": {"evalue_cutoff": evalue, "identity_cutoff": 0, "sequence_type": "protein", "value": seq}},
         "return_type": "polymer_entity",
         "request_options": {"scoring_strategy": "sequence", "paginate": {"start": 0, "rows": 10},
                             "results_content_type": ["experimental"]}}
    req = urllib.request.Request(URL, data=json.dumps(q).encode(), headers={'Content-Type': 'application/json'})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if r.status == 204:
                    return 0, []
                d = json.loads(r.read() or b'{}')
                hits = []
                for h in d.get('result_set', []):
                    ctx = (h.get('services') or [{}])[0].get('nodes', [{}])[0].get('match_context', [{}])[0]
                    hits.append((h['identifier'], ctx.get('evalue'), ctx.get('sequence_identity')))
                return d.get('total_count', 0), hits
        except Exception as e:
            time.sleep(3 * (attempt + 1)); err = e
    raise err


if __name__ == '__main__':
    with open(sys.argv[2], 'w') as out:
        for line in open(sys.argv[1]):
            acc, seq = line.split()[:2]
            n, hits = search(seq)
            best = hits[0] if hits else ('-', '-', '-')
            out.write('%s\t%d\t%s\t%s\t%s\n' % (acc, n, *best)); out.flush()
            print(acc, n, best, flush=True)
