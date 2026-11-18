"""Divergence-site test: at sites where a post-cutoff protein differs from its closest pre-cutoff relative,
does Fable output the true (new) residue or the old relative's residue?
Baselines from the profile of all other old relatives at that column: f_true, f_rel (what a model that knows the
family but not this particular relative would be expected to give, if it sampled from the family profile)."""
import json, sys, csv, collections, random
rows = [json.loads(l) for l in open(sys.argv[1])]
rows = [r for r in rows if 'error' not in r and r.get('valid')]
def ci(k, n):
    import math
    if n == 0: return (0, 0)
    p = k / n; z = 1.96; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (c - h), 100 * (c + h))
out = []
for site in ['div', 'cons']:
    rs = [r for r in rows if r['site'] == site]
    if not rs: continue
    n = len(rs)
    k_true = sum(r['pred'][0] == r['true_res'] for r in rs)
    k_rel = sum(r['pred'][0] == r['rel_res'] for r in rs)
    k_maj = sum(r['pred'][0] == r['majority'] for r in rs if r['majority'])
    ft = sum(float(r['f_true'] or 0) for r in rs) / n; fr = sum(float(r['f_rel'] or 0) for r in rs) / n
    maj_true = sum(r['majority'] == r['true_res'] for r in rs) / n; maj_rel = sum(r['majority'] == r['rel_res'] for r in rs) / n
    print(f'== {site} sites: n={n} proteins={len({r["accession"] for r in rs})}')
    print('  Fable says true (new) residue : %5.1f%%  [%4.1f, %4.1f]' % (100 * k_true / n, *ci(k_true, n)))
    if site == 'div':
        print('  Fable says old relative residue: %5.1f%%  [%4.1f, %4.1f]' % (100 * k_rel / n, *ci(k_rel, n)))
        print('  Fable says neither             : %5.1f%%' % (100 * (n - k_true - k_rel) / n))
    print('  profile of other old relatives: mean freq of true residue %.1f%%, of relative residue %.1f%%' % (100 * ft, 100 * fr))
    print('  majority of other old relatives = true %.1f%%, = relative %.1f%%' % (100 * maj_true, 100 * maj_rel))
    print('  Fable says the majority residue: %.1f%%' % (100 * k_maj / n))
    # split by whether the relative's residue is the family majority
    if site == 'div':
        for lab, sel in [('relative residue = majority of others', lambda r: r['majority'] == r['rel_res']),
                         ('relative residue != majority of others', lambda r: r['majority'] != r['rel_res']),
                         ('relative residue rare in others (<20%)', lambda r: r['f_rel'] != '' and float(r['f_rel']) < 0.2)]:
            ss = [r for r in rs if sel(r)]
            if not ss: continue
            m = len(ss)
            print('   %-42s n=%3d  Fable=true %5.1f%%  Fable=relative %5.1f%%  | profile f_true %.1f%% f_rel %.1f%%' % (
                lab, m, 100 * sum(r['pred'][0] == r['true_res'] for r in ss) / m, 100 * sum(r['pred'][0] == r['rel_res'] for r in ss) / m,
                100 * sum(float(r['f_true'] or 0) for r in ss) / m, 100 * sum(float(r['f_rel'] or 0) for r in ss) / m))
# continuation: does Fable's 10-residue reply match the relative's continuation better than the new protein's?
