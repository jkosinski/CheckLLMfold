"""Contact prediction and modelling in one reply, without reasoning: no tools, no thinking (effort low,
CLAUDE_CODE_DISABLE_THINKING=1; thinking tokens recorded per call), and the reply must be the answer only.
Same headless setup as scripts/recognize.py.

Usage: oneshot.py MODE OUT.jsonl [--reps 3] [--jobs 6] [--effort low]
MODE: contacts (ranked list of pairs with |i-j| >= 6, as many as in the main protocol) or coords (C-alpha coordinates).
Proteins: data/recognition/contact7.tsv; number of pairs from data/targets.tsv and data/famous_targets.tsv."""
import os, re, json, csv, subprocess, tempfile, argparse, concurrent.futures as cf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETTINGS = os.path.join(ROOT, 'runs', 'sandbox_settings.json')


def fasta(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def prompt(mode, seq, npairs):
    L = len(seq)
    if mode == 'contacts':
        return (f"Below is a protein sequence in one-letter amino acid codes ({L} residues, numbered from 1). Predict "
                "which residue pairs are in contact in its folded structure: C-beta atoms (C-alpha for glycine) within "
                f"8 A, at sequence separation |i - j| >= 6. Reply with exactly {npairs} pairs, ranked from most to least "
                "confident, one pair per line as two numbers 'i j' with i < j, and nothing else. Do not think it over "
                "or work it out first: write the answer straight away.\n\n"
                f">sequence\n{fasta(seq)}")
    return (f"Below is a protein sequence in one-letter amino acid codes ({L} residues, numbered from 1). Predict its "
            "folded 3D structure. Reply with the C-alpha coordinates in Angstrom for all residues in order, one line per "
            f"residue as 'number x y z' (lines 1 to {L}), and nothing else. Do not think it over or work it out "
            f"first: write the answer straight away.\n\n>sequence\n{fasta(seq)}")


def parse(mode, text, L):
    out = []
    for line in (text or '').splitlines():
        nums = re.findall(r'-?\d+(?:\.\d+)?', line)
        if mode == 'contacts' and len(nums) == 2 and '.' not in line:
            i, j = sorted(int(x) for x in nums)
            if 1 <= i and j <= L and j - i >= 6:
                out.append([i, j])
        elif mode == 'coords' and len(nums) == 4:
            out.append([int(float(nums[0]))] + [float(x) for x in nums[1:]])
    return out


def ask(q, mode, model, effort):
    env = dict(os.environ, CLAUDE_CODE_DISABLE_THINKING='1', MAX_THINKING_TOKENS='0')  # Fable 5.2 thinks anyway
    with tempfile.TemporaryDirectory(prefix='q_') as work:
        p = subprocess.run(
            ['claude', '-p', '--model', model, '--effort', effort, '--tools', '', '--strict-mcp-config',
             '--setting-sources', 'project', '--settings', SETTINGS, '--disable-slash-commands',
             '--system-prompt', 'You are a helpful assistant.', '--no-session-persistence',
             '--output-format', 'json'],
            input=prompt(mode, q['sequence'], q['npairs']), capture_output=True, text=True, cwd=work, timeout=1800, env=env)
    try:
        o = json.loads(p.stdout)
    except Exception:
        return dict(q, mode=mode, error=(p.stdout + p.stderr)[-500:])
    text = o.get('result') or ''
    u = o.get('usage') or {}
    return dict(q, mode=mode, model=model, effort=effort, text=text, parsed=parse(mode, text, len(q['sequence'])),
                thinking_tokens=(u.get('output_tokens_details') or {}).get('thinking_tokens'),
                output_tokens=u.get('output_tokens'), is_error=o.get('is_error'), cost_usd=o.get('total_cost_usd'),
                duration_s=(o.get('duration_ms') or 0) / 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=('contacts', 'coords')); ap.add_argument('out')
    ap.add_argument('--model', default='claude-fable-5-2'); ap.add_argument('--effort', default='low')
    ap.add_argument('--reps', type=int, default=3); ap.add_argument('--jobs', type=int, default=6)
    ap.add_argument('--only', default='')
    a = ap.parse_args()
    npairs = {}
    for f in ('targets.tsv', 'famous_targets.tsv'):
        for r in csv.reader(open(os.path.join(ROOT, 'data', f)), delimiter='\t'):
            npairs[r[3]] = int(r[2])
    npairs['MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'] = 286
    qs = []
    for q in csv.DictReader(open(os.path.join(ROOT, 'data/recognition/contact7.tsv')), delimiter='\t'):
        if a.only and a.only not in q['qid']:
            continue
        for r in range(1, a.reps + 1):
            qs.append(dict(qid='%s|%s|rep%d' % (q['qid'].split('|')[1], a.mode, r), target=q['qid'].split('|')[1],
                           set=q['set'], sequence=q['sequence'], npairs=npairs[q['sequence']]))
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            o = json.loads(l)
            if 'error' not in o and not o.get('is_error'):
                done.add(o['qid'])
    qs = [q for q in qs if q['qid'] not in done]
    print('to ask:', len(qs), 'already done:', len(done), flush=True)
    with cf.ThreadPoolExecutor(a.jobs) as ex, open(a.out, 'a') as f:
        for r in ex.map(lambda q: ask(q, a.mode, a.model, a.effort), qs):
            f.write(json.dumps(r) + '\n'); f.flush()
            print(r['qid'], 'ERROR ' + (r.get('error') or r.get('text', ''))[-200:] if 'error' in r or r.get('is_error')
                  else 'parsed %d thinking %s out %s' % (len(r['parsed']), r['thinking_tokens'], r['output_tokens']), flush=True)


if __name__ == '__main__':
    main()
