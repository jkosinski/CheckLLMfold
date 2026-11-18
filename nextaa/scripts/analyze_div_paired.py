"""Pair each divergence-test site (new protein, prefix up to k) with the same site on its closest old relative
(relative's own prefix up to the aligned position). Usage: analyze_div_paired.py DIVDIR NEW.jsonl REL.jsonl"""
import sys, json, csv, collections, math
d, fnew, frel = sys.argv[1:4]
q = {r['qid']: r for r in csv.DictReader(open(f'{d}/q_div.tsv'), delimiter='\t')}
new = {o['qid']: o for o in map(json.loads, open(fnew)) if 'error' not in o and o.get('valid')}
rel = {}
for o in map(json.loads, open(frel)):
    if 'error' in o or not o.get('valid'): continue
    rel[(o['new_accession'], o['site'], o['condition'], o['accession'], int(o['cut']))] = o
relq = list(csv.DictReader(open(f'{d}/q_rel.tsv'), delimiter='\t'))
# map new-site qid -> relative-site question (same order of construction)
pairs = []
qs = list(csv.DictReader(open(f'{d}/q_div.tsv'), delimiter='\t'))
rq_real = [r for r in relq if r['condition'] == 'real']
assert len(rq_real) == len(qs)
for a, b in zip(qs, rq_real):
    assert a['accession'] == b['new_accession'] and a['site'] == b['site']
    n = new.get(a['qid']); r = rel.get((b['new_accession'], b['site'], 'real', b['accession'], int(b['cut'])))
    if n and r: pairs.append((a, n, r))
def pct(k, n): return '%5.1f%% (%d/%d)' % (100 * k / n if n else 0, k, n)
print('pairs', len(pairs))
# relative memorisation overall
R = [r for a, n, r in pairs]
print('Fable on the OLD relative itself: first residue right', pct(sum(r['first_correct'] for r in R), len(R)),
      '; 10/10 exact', sum(r['identity'] == 1 for r in R), '; >=5 first right', sum(r['pred'][:5] == r['target'][:5] for r in R))
shuf = [o for k, o in rel.items() if k[2] == 'shuffled']
if shuf: print('   shuffled control on relatives:', pct(sum(o['first_correct'] for o in shuf), len(shuf)))
for site in ['div', 'cons']:
    P = [(a, n, r) for a, n, r in pairs if a['site'] == site]
    print(f'\n== {site} sites (new residue {"!=" if site == "div" else "=="} relative residue), n={len(P)}')
    for lab, sel in [('relative recalled (Fable right on the relative here)', lambda r: r['first_correct']),
                     ('relative not recalled', lambda r: not r['first_correct'])]:
        S = [(a, n, r) for a, n, r in P if sel(r)]
        if not S: continue
        say_true = sum(n['pred'][0] == a['true_res'] for a, n, r in S)
        say_rel = sum(n['pred'][0] == a['rel_res'] for a, n, r in S)
        ft = sum(float(a['f_true'] or 0) for a, n, r in S) / len(S); fr = sum(float(a['f_rel'] or 0) for a, n, r in S) / len(S)
        if site == 'div':
            print(f'  {lab:52s} n={len(S):3d}  on new protein Fable says: new residue {pct(say_true, len(S))}, relative residue {pct(say_rel, len(S))} | family profile: new {100*ft:.0f}%, relative {100*fr:.0f}%')
        else:
            print(f'  {lab:52s} n={len(S):3d}  on new protein Fable right {pct(say_true, len(S))} | family profile {100*ft:.0f}%')
# by relative creation era
print()
for lab, sel in [('relative in UniProtKB before 2025', lambda a: a['rel_created'] < '2025/01/01'), ('relative added 2025-2026', lambda a: a['rel_created'] >= '2025/01/01')]:
    S = [r for a, n, r in pairs if sel(a)]
    if S: print(f'{lab:32s} Fable right on relative itself {pct(sum(r["first_correct"] for r in S), len(S))}')
