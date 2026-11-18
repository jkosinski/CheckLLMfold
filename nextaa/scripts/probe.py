"""Next-amino-acid probe, replicating Tu & Chen's chat prompt (LLM-bioseq-probe, bioseq_probe.probe.prompt) verbatim.
Runs Claude through headless Claude Code (`claude -p`): no tools, no MCP, CLAUDE.md files excluded, empty working dir.

Usage: probe.py QUESTIONS.tsv OUT.jsonl [--model claude-fable-5-2] [--effort low] [--jobs 6]
QUESTIONS.tsv columns: qid, prefix, target (10 residues), plus any metadata columns (copied to the output).
Questions already in OUT.jsonl are skipped, so a run can be resumed."""
import sys, os, json, csv, re, subprocess, tempfile, argparse, concurrent.futures as cf, time

K = 10
STANDARD = 'ACDEFGHIKLMNPQRSTVWY'
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETTINGS = os.path.join(os.path.dirname(ROOT), 'runs', 'sandbox_settings.json')


def fasta(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def prompt(prefix):
    return ("Below is the beginning of a protein sequence in one-letter amino acid codes. "
            f"Predict how it continues. Reply with exactly the next {K} amino acids as one-letter "
            "codes on one line, and nothing else.\n\n"
            f">sequence\n{fasta(prefix)}")


def parse(text, prefix):
    """Last line that is only standard residues (spaces ignored); a restatement of the prefix's end is removed."""
    seqs = []
    for line in (text or '').splitlines():
        s = re.sub(r'[\s`*]', '', line)
        if s and re.fullmatch(f'[{STANDARD}]+', s):
            seqs.append(s)
    if not seqs:
        return ''
    s = seqs[-1]
    for n in range(min(len(s), len(prefix)), 14, -1):  # reply opening with >=15 residues of the prefix's end
        if prefix.endswith(s[:n]):
            return s[n:]
    return s


def ask(q, model, effort):
    with tempfile.TemporaryDirectory(prefix='q_') as work:
        p = subprocess.run(
            ['claude', '-p', '--model', model, '--effort', effort, '--tools', '', '--strict-mcp-config',
             '--setting-sources', 'project', '--settings', SETTINGS, '--disable-slash-commands',
             '--system-prompt', 'You are a helpful assistant.', '--no-session-persistence',
             '--output-format', 'json'],
            input=prompt(q['prefix']), capture_output=True, text=True, cwd=work, timeout=600)
    try:
        o = json.loads(p.stdout)
    except Exception:
        return dict(q, error=(p.stdout + p.stderr)[-500:])
    text = o.get('result') or ''
    pred = parse(text, q['prefix'])
    t = q['target']
    return dict(q, model=model, effort=effort, text=text, pred=pred, valid=bool(pred),
                first_correct=bool(pred) and pred[0] == t[0],
                identity=sum(a == b for a, b in zip(pred, t)) / len(t),
                is_error=o.get('is_error'), cost_usd=o.get('total_cost_usd'), usage=o.get('usage'),
                duration_s=(o.get('duration_ms') or 0) / 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('questions'); ap.add_argument('out')
    ap.add_argument('--model', default='claude-fable-5-2'); ap.add_argument('--effort', default='low')
    ap.add_argument('--jobs', type=int, default=6); ap.add_argument('--limit', type=int, default=0)
    a = ap.parse_args()
    qs = list(csv.DictReader(open(a.questions), delimiter='\t'))
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            o = json.loads(l)
            if 'error' not in o and o.get('model') == a.model: done.add(o['qid'])
    qs = [q for q in qs if q['qid'] not in done]
    if a.limit: qs = qs[:a.limit]
    print('to ask:', len(qs), 'already done:', len(done), flush=True)
    n = 0
    with cf.ThreadPoolExecutor(a.jobs) as ex, open(a.out, 'a') as f:
        for r in ex.map(lambda q: ask(q, a.model, a.effort), qs):
            f.write(json.dumps(r) + '\n'); f.flush(); n += 1
            if 'error' in r: print('ERROR', r['qid'], r['error'][-200:], flush=True)
            elif n % 25 == 0: print('done', n, 'of', len(qs), flush=True)
    print('finished', n, flush=True)


if __name__ == '__main__':
    main()
