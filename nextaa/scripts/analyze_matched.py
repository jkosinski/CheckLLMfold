"""Gap (real - shuffled first-residue accuracy) for our Fable 5.2 runs: well-known, pre-cutoff UniRef50 singletons,
post-cutoff D set; mean over families of per-family gap with 95% bootstrap CI, plus pooled accuracies.
Also the same subsets from Tu & Chen's published Fable 5.2 replies, for comparison."""
import json, gzip, csv, collections, random, sys

def fam_gap(rows, key='pfam'):
    fam = collections.defaultdict(lambda: {'real': [], 'shuffled': []})
    for o in rows: fam[o[key] or o['accession']][o['condition']].append(bool(o['first_correct']))
    fam = {k: v for k, v in fam.items() if v['real'] and v['shuffled']}
    g = {k: sum(v['real']) / len(v['real']) - sum(v['shuffled']) / len(v['shuffled']) for k, v in fam.items()}
    keys = sorted(g); m = sum(g.values()) / len(g)
    rnd = random.Random(0); bs = sorted(sum(g[rnd.choice(keys)] for _ in keys) / len(keys) for _ in range(4000))
    R = [o['first_correct'] for o in rows if o['condition'] == 'real']; S = [o['first_correct'] for o in rows if o['condition'] == 'shuffled']
    return dict(families=len(keys), proteins=len({o['accession'] for o in rows}), real_n=len(R), shuf_n=len(S),
                real=100 * sum(R) / len(R), shuffled=100 * sum(S) / len(S), gap=100 * m, lo=100 * bs[100], hi=100 * bs[3899])

ours = [json.loads(l) for f in ['results/famous_fable52.jsonl', 'results/matched_fable52.jsonl'] for l in open(f)]
ours = [o for o in ours if 'error' not in o and o.get('valid')]
ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open('data/main_uniref50.tsv'), delimiter='\t')}
qids = {r['qid'] for f in ['data/q_pre_single.tsv', 'data/q_post_D.tsv'] for r in csv.DictReader(open(f), delimiter='\t')}
tu = []
for f, s in [('data/tu_outputs/confirm/claude-fable-5-2.jsonl.gz', 'pre_single'), ('data/tu_outputs/d_set/claude-fable-5-2.jsonl.gz', 'post_D')]:
    for l in gzip.open(f, 'rt'):
        o = json.loads(l)
        if o.get('valid') is False: continue
        if s == 'pre_single' and ur50[o['accession']] != 1: continue
        o['set'] = s; tu.append(o)
out = csv.writer(open('results/matched_summary.tsv', 'w'), delimiter='\t')
out.writerow(['source', 'set', 'families', 'proteins', 'real_n', 'shuf_n', 'real_acc', 'shuffled_acc', 'gap', 'ci_lo', 'ci_hi'])
for src, rows in [('this run (claude -p, effort low)', ours), ('Tu & Chen replies (OpenRouter)', tu)]:
    for s in ['famous', 'pre_single', 'post_D']:
        rs = [o for o in rows if o['set'] == s]
        if not rs: continue
        r = fam_gap(rs, 'pfam' if s != 'famous' else 'accession')
        out.writerow([src, s, r['families'], r['proteins'], r['real_n'], r['shuf_n'], '%.1f' % r['real'], '%.1f' % r['shuffled'], '%.1f' % r['gap'], '%.1f' % r['lo'], '%.1f' % r['hi']])
        print('%-34s %-11s fam %3d prot %3d  n %3d/%3d  real %5.1f  shuf %4.1f  gap %5.1f [%5.1f, %5.1f]' % (src, s, r['families'], r['proteins'], r['real_n'], r['shuf_n'], r['real'], r['shuffled'], r['gap'], r['lo'], r['hi']))
# exact 10-residue recall
for s in ['famous', 'pre_single', 'post_D']:
    rs = [o for o in ours if o['set'] == s and o['condition'] == 'real']
    if rs: print('exact 10/10 continuation', s, sum(o['identity'] == 1 for o in rs), 'of', len(rs))
print('cost so far $%.2f' % sum(o.get('cost_usd') or 0 for o in ours))
