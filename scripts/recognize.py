"""Recognition probe: give Fable a full protein sequence (no tools) and ask what protein it is.
Same headless setup as nextaa/scripts/probe.py: no tools, no MCP, CLAUDE.md files excluded, empty working dir.

Usage: recognize.py QUESTIONS.tsv OUT.jsonl [--model claude-fable-5-2] [--effort low] [--reps 1] [--jobs 6]
QUESTIONS.tsv columns: qid, sequence, plus any metadata columns (copied to the output).
Questions already in OUT.jsonl are skipped, so a run can be resumed."""
import os, re, json, csv, subprocess, tempfile, argparse, concurrent.futures as cf

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SETTINGS = os.path.join(ROOT, 'runs', 'sandbox_settings.json')
FIELDS = ('name', 'accession', 'organism', 'pfam', 'fold', 'pdb')


def fasta(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def prompt(seq):
    return ("Below is a protein sequence in one-letter amino acid codes. Identify it as precisely as you can. "
            "Answer in exactly this format, one field per line, and write unknown for any field you cannot give:\n"
            "name: <protein name>\n"
            "accession: <UniProt or RefSeq accession>\n"
            "organism: <source organism>\n"
            "pfam: <Pfam accession and name>\n"
            "fold: <fold or structural class>\n"
            "pdb: <PDB IDs of this protein or a close relative>\n\n"
            f">sequence\n{fasta(seq)}")


def parse(text):
    out = {}
    for line in (text or '').splitlines():
        m = re.match(r'\s*[*_]*(name|accession|organism|pfam|fold|pdb)[*_]*\s*:\s*(.*)', line, re.I)
        if m and m.group(1).lower() not in out:
            out[m.group(1).lower()] = m.group(2).strip().strip('*_ ')
    return {k: out.get(k, '') for k in FIELDS}


def ask(q, model, effort):
    with tempfile.TemporaryDirectory(prefix='q_') as work:
        p = subprocess.run(
            ['claude', '-p', '--model', model, '--effort', effort, '--tools', '', '--strict-mcp-config',
             '--setting-sources', 'project', '--settings', SETTINGS, '--disable-slash-commands',
             '--system-prompt', 'You are a helpful assistant.', '--no-session-persistence',
             '--output-format', 'json'],
            input=prompt(q['sequence']), capture_output=True, text=True, cwd=work, timeout=900)
    try:
        o = json.loads(p.stdout)
    except Exception:
        return dict(q, error=(p.stdout + p.stderr)[-500:])
    text = o.get('result') or ''
    return dict(q, model=model, effort=effort, text=text, answer=parse(text), is_error=o.get('is_error'),
                cost_usd=o.get('total_cost_usd'), usage=o.get('usage'), duration_s=(o.get('duration_ms') or 0) / 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('questions'); ap.add_argument('out')
    ap.add_argument('--model', default='claude-fable-5-2'); ap.add_argument('--effort', default='low')
    ap.add_argument('--reps', type=int, default=1); ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args()
    qs = []
    for q in csv.DictReader(open(a.questions), delimiter='\t'):
        for r in range(1, a.reps + 1):
            qs.append(dict(q, qid='%s|rep%d' % (q['qid'], r)) if a.reps > 1 else q)
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            o = json.loads(l)
            if 'error' not in o and not o.get('is_error') and o.get('model') == a.model and o.get('effort') == a.effort:
                done.add(o['qid'])
    qs = [q for q in qs if q['qid'] not in done]
    print('to ask:', len(qs), 'already done:', len(done), flush=True)
    n = 0
    with cf.ThreadPoolExecutor(a.jobs) as ex, open(a.out, 'a') as f:
        for r in ex.map(lambda q: ask(q, a.model, a.effort), qs):
            f.write(json.dumps(r) + '\n'); f.flush(); n += 1
            if 'error' in r or r.get('is_error'):
                print('ERROR', r['qid'], (r.get('error') or r.get('text', ''))[-200:], flush=True)
            elif n % 20 == 0:
                print('done', n, 'of', len(qs), flush=True)
    print('finished', n, flush=True)


if __name__ == '__main__':
    main()
