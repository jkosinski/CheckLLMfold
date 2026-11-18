"""Re-analyse Tu & Chen's published replies: next-residue gap (real - shuffled) on the main set split by
UniRef50 cluster size, against the post-cutoff D set. Gap = mean over families of per-family gap (their statistic),
with a 95% bootstrap CI over families; also pooled per-cut accuracy."""
import gzip, json, csv, random, collections, sys
OUT = 'data/tu_outputs'
ur50 = {r['accession']: int(r['member_count']) for r in csv.DictReader(open('data/main_uniref50.tsv'), delimiter='\t')}
def stratum(a):
    c = ur50.get(a)
    if c is None: return 'post-cutoff (D set)'
    return 'UniRef50 singleton' if c == 1 else ('UniRef50 2-5' if c <= 5 else ('UniRef50 6-50' if c <= 50 else 'UniRef50 >50'))
def load(path):
    rows = []
    for l in gzip.open(path, 'rt'):
        o = json.loads(l)
        if o.get('valid') is False: continue
        rows.append(o)
    return rows
def gap(rows, families=None):
    fam = collections.defaultdict(lambda: {'real': [], 'shuffled': []})
    for o in rows: fam[o['pfam']][o['condition']].append(bool(o['first_correct']))
    fam = {k: v for k, v in fam.items() if v['real'] and v['shuffled']}
    g = {k: sum(v['real']) / len(v['real']) - sum(v['shuffled']) / len(v['shuffled']) for k, v in fam.items()}
    keys = sorted(g); m = sum(g[k] for k in keys) / len(keys)
    rnd = random.Random(0); bs = []
    for _ in range(4000):
        s = [g[rnd.choice(keys)] for _ in keys]; bs.append(sum(s) / len(s))
    bs.sort()
    real = [x for o in rows if o['condition'] == 'real' for x in [o['first_correct']]]
    shuf = [x for o in rows if o['condition'] == 'shuffled' for x in [o['first_correct']]]
    return dict(fams=len(keys), prots=len({o['accession'] for o in rows}), cuts=len(real),
                gap=100 * m, lo=100 * bs[100], hi=100 * bs[3899],
                real=100 * sum(real) / len(real), shuf=100 * sum(shuf) / len(shuf))
models = [('Claude Fable 5.2', 'confirm/claude-fable-5-2', 'd_set/claude-fable-5-2'),
          ('Claude Opus 5.5', 'confirm/claude-opus-5-5', 'd_set/claude-opus-5-5'),
          ('ProGen3-3B', 'plm_confirm/progen3-3b__prefix', 'plm_d/progen3-3b__prefix'),
          ('ProtGPT3-10B', 'plm_confirm/protgpt3-10B__prefix', 'plm_d/protgpt3-10B__prefix')]
order = ['UniRef50 >50', 'UniRef50 6-50', 'UniRef50 2-5', 'UniRef50 singleton', 'post-cutoff (D set)']
w = csv.writer(open('results/tu_reanalysis_by_uniref50.tsv', 'w'), delimiter='\t')
w.writerow(['model', 'stratum', 'families', 'proteins', 'cuts', 'gap', 'ci_lo', 'ci_hi', 'real_acc', 'shuffled_acc'])
for name, main, d in models:
    rows = load(f'{OUT}/{main}.jsonl.gz') + load(f'{OUT}/{d}.jsonl.gz')
    by = collections.defaultdict(list)
    for o in rows: by[stratum(o['accession'])].append(o)
    # families shared by D set and the main-set singletons
    print('==', name)
    for s in order:
        r = gap(by[s])
        w.writerow([name, s, r['fams'], r['prots'], r['cuts'], '%.1f' % r['gap'], '%.1f' % r['lo'], '%.1f' % r['hi'], '%.1f' % r['real'], '%.1f' % r['shuf']])
        print('%-22s fam %3d prot %3d cuts %4d  gap %5.1f [%5.1f, %5.1f]  real %4.1f shuf %4.1f' % (s, r['fams'], r['prots'], r['cuts'], r['gap'], r['lo'], r['hi'], r['real'], r['shuf']))
