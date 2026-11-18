"""Score the recognition probe (scripts/recognize.py) -> results/recognition/summary.tsv.
Accession: the true accession appears in the accession field. Pfam: the true Pfam accession or short name (e.g. DUF222)
appears in the pfam or name field. Genus: the true genus appears in the organism field. Chance level: the same
scores with the answers shuffled among the proteins of a set (mean of 1000 shuffles), which credits a model that
gives the same common family or genus to many proteins."""
import os, re, csv, json, random, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results', 'recognition')
short = {r['pfam']: r['short'] for r in csv.DictReader(open(os.path.join(ROOT, 'data/recognition/pfam_names.tsv')), delimiter='\t')}


def load(p):
    rows = {}
    for l in open(p):
        o = json.loads(l)
        if 'error' not in o and not o.get('is_error'):
            rows[o['qid']] = o
    return list(rows.values())


def word(s, text):
    return bool(s) and re.search(r'(?<![A-Za-z0-9])' + re.escape(s) + r'(?![A-Za-z0-9])', text or '', re.I) is not None


def genus(org):
    w = [x for x in re.split(r'\s+', org or '') if x and x.lower() != 'candidatus']
    return w[0].strip('[]') if w else ''


def m_acc(a, t): return word(t['accession'], a['accession'])
def m_pfam(a, t):
    if not t.get('pfam'):
        d = re.search(r'DUF\d+', t.get('truth', ''))
        return bool(d) and (word(d.group(0), a['pfam']) or word(d.group(0), a['name']))
    return word(t['pfam'], a['pfam']) or word(short.get(t['pfam'], ''), a['pfam']) or word(short.get(t['pfam'], ''), a['name'])
def m_genus(a, t): return word(genus(t['organism']), a['organism'])
def m_pdb(a, t): return word(t.get('pdb', ''), a['pdb'])
def unknown(v): return not v or v.strip().lower().startswith('unknown')


def chance(rows, f, n=1000, seed=1):
    rng = random.Random(seed); ans = [r['answer'] for r in rows]; tot = 0
    for _ in range(n):
        rng.shuffle(ans); tot += sum(f(a, t) for a, t in zip(ans, rows))
    return tot / n / len(rows)


out = []
c7 = load(os.path.join(R, 'contact7_low.jsonl')) + load(os.path.join(R, 'contact7_high.jsonl'))
print('Contact-test proteins (3 repeats per effort)')
by = collections.defaultdict(list)
for o in c7:
    by[(o['qid'].rsplit('|', 1)[0], o['effort'])].append(o)
for (q, eff), os_ in sorted(by.items()):
    t = os_[0]
    row = dict(set='contact7', group=q.split('|')[1], effort=eff, n=len(os_),
               accession=sum(m_acc(o['answer'], t) for o in os_), pfam=sum(m_pfam(o['answer'], t) for o in os_),
               genus=sum(m_genus(o['answer'], t) for o in os_), pdb=sum(m_pdb(o['answer'], t) for o in os_) if t['pdb'] else '',
               name_unknown=sum(unknown(o['answer']['name']) for o in os_),
               example_name=os_[0]['answer']['name'][:70])
    out.append(row)
    print('  %-12s %-4s acc %d/%d pfam %d genus %d pdb %s | %s' % (row['group'], eff, row['accession'], row['n'], row['pfam'],
                                                               row['genus'], row['pdb'], row['example_name']))

na = load(os.path.join(R, 'nextaa_low.jsonl'))
print('\nNext-residue proteins (effort low, one question each)')
print('  %-10s %3s %-14s %-14s %-14s %s' % ('set', 'n', 'accession', 'pfam (chance)', 'genus (chance)', 'unknown: acc/pfam/organism'))
for s in ('famous', 'pre_big', 'pre_single', 'post_D'):
    rows = [o for o in na if o['set'] == s]
    if not rows:
        continue
    f = lambda fn: sum(fn(o['answer'], o) for o in rows) / len(rows)
    row = dict(set='nextaa', group=s, effort='low', n=len(rows), accession=round(f(m_acc), 3),
               pfam=round(f(m_pfam), 3) if s != 'famous' else '', pfam_chance=round(chance(rows, m_pfam), 3) if s != 'famous' else '',
               genus=round(f(m_genus), 3) if s != 'famous' else '', genus_chance=round(chance(rows, m_genus), 3) if s != 'famous' else '',
               acc_unknown=round(sum(unknown(o['answer']['accession']) for o in rows) / len(rows), 3),
               pfam_unknown=round(sum(unknown(o['answer']['pfam']) for o in rows) / len(rows), 3),
               organism_unknown=round(sum(unknown(o['answer']['organism']) for o in rows) / len(rows), 3))
    out.append(row)
    print('  %-10s %3d %-14s %-14s %-14s %s/%s/%s' % (s, row['n'], row['accession'], '%s (%s)' % (row['pfam'], row['pfam_chance']),
                                                    '%s (%s)' % (row['genus'], row['genus_chance']), row['acc_unknown'],
                                                    row['pfam_unknown'], row['organism_unknown']))
cols = ['set', 'group', 'effort', 'n', 'accession', 'pfam', 'pfam_chance', 'genus', 'genus_chance', 'pdb',
        'name_unknown', 'acc_unknown', 'pfam_unknown', 'organism_unknown', 'example_name']
with open(os.path.join(R, 'summary.tsv'), 'w') as fh:
    w = csv.DictWriter(fh, cols, delimiter='\t', restval=''); w.writeheader(); w.writerows(out)
