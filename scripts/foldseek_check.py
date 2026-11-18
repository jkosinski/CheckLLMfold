"""Submit structures to the public Foldseek server (search.foldseek.com) against PDB100 and save hits.
Usage: foldseek_check.py OUTDIR MODE(3diaa|tmalign) file1.pdb [file2.pdb ...]"""
import sys, os, json, time, subprocess

API = 'https://search.foldseek.com/api'


def curl(args):
    r = subprocess.run(['curl', '-sS', '-m', '120'] + args, capture_output=True, text=True)
    return r.stdout


def submit(path, mode):
    out = curl(['-X', 'POST', '-F', 'q=@' + path, '--form-string', 'mode=' + mode,
                '--form-string', 'database[]=pdb100', API + '/ticket'])
    return json.loads(out)


def wait(tid):
    for _ in range(120):
        st = json.loads(curl([API + '/ticket/' + tid]))
        if st.get('status') == 'COMPLETE':
            return True
        if st.get('status') == 'ERROR':
            return False
        time.sleep(5)
    return False


def main():
    outdir, mode, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    os.makedirs(outdir, exist_ok=True)
    for f in files:
        acc = os.path.basename(f).split('.')[0]
        dest = os.path.join(outdir, '%s.%s.json' % (acc, mode))
        if os.path.exists(dest):
            continue
        t = submit(f, mode)
        tid = t.get('id')
        if not tid or not wait(tid):
            print(acc, 'FAILED', t, flush=True); continue
        res = json.loads(curl([API + '/result/' + tid + '/0']))
        json.dump(res, open(dest, 'w'))
        hits = []
        for db in res.get('results', []):
            for aln in (db.get('alignments') or []):
                for a in (aln if isinstance(aln, list) else [aln]):
                    hits.append(a)
        hits.sort(key=lambda a: -float(a.get('prob', 0) or 0))
        top = hits[0] if hits else {}
        print('%s %s nhits=%d top=%s prob=%s eval=%s tm=%s' % (acc, mode, len(hits), top.get('target', '-')[:40],
              top.get('prob'), top.get('eval'), top.get('alntmscore')), flush=True)
        time.sleep(3)


if __name__ == '__main__':
    main()
