"""Score all runs listed in runs/jobs.tsv (plus Sergey's shared chat) and write results/summary.tsv."""
import os, sys, json, subprocess
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = {'ubiquitin_1UBQ': ('data/1UBQ/1ubq.pdb', 'A')}
for line in open(os.path.join(ROOT, 'data/targets.tsv')):
    a = line.split('\t')[0]
    REF[a] = ('data/%s/%s.pdb' % (a, a), None)
for line in open(os.path.join(ROOT, 'data/famous.tsv')):
    name, f, ch = line.split()
    REF[name] = (f, ch)

rows = []
runs = [l.rstrip('\n').split('\t') for l in open(os.path.join(ROOT, 'runs/jobs.tsv')) if l.strip()]
runs.append(['sergey_shared_chat', 'ubiquitin_1UBQ', 'claude.ai', 'claude.ai'])
for job, target, rep, effort in runs:
    d = os.path.join(ROOT, 'runs/jobs', job) if job.startswith('job_') else os.path.join(ROOT, 'runs', job)
    if not os.path.exists(os.path.join(d, 'final.txt')):
        print('skip (not finished)', job, target); continue
    ref, ch = REF[target]
    args = [sys.executable, '-I', os.path.join(ROOT, 'scripts/evaluate.py'), d, os.path.join(ROOT, ref)] + ([ch] if ch else [])
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode:
        print('FAILED', job, target, r.stderr[-500:]); continue
    m = json.load(open(os.path.join(d, 'metrics.json')))
    cost = dur = None
    tj = os.path.join(d, 'transcript.jsonl')
    if os.path.exists(tj):
        for l in open(tj):
            try:
                o = json.loads(l)
            except Exception:
                continue
            if o.get('type') == 'result':
                cost, dur = o.get('total_cost_usd'), (o.get('duration_ms') or 0) / 1000
    rows.append(dict(job=job, target=target, rep=rep, cost_usd=cost, duration_s=dur, **m))

os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
keys = ['job', 'target', 'rep', 'L', 'rounds', 'n_final_contacts', 'native_contacts_sep6', 'chance_precision_sep6',
        'P_sep6_topL5', 'P_sep6_topL2', 'P_sep6_topL', 'P_sep6_top2L', 'P_sep24_topL5', 'P_sep24_topL2', 'P_sep24_topL',
        'final_P_sep6', 'final_P_sep24', 'final_recall_sep6',
        'final_P_local3_5', 'final_R_local3_5', 'final_P_short6_11', 'final_R_short6_11',
        'final_P_medium12_23', 'final_R_medium12_23', 'final_P_long24', 'final_R_long24',
        'TM_model_vs_ref', 'TM_mirror_vs_ref',
        'TM_ceiling_native_contacts', 'TM_ceiling_mirror', 'TM_best_hand', 'TM_ceiling_best_hand', 'cost_usd', 'duration_s']
with open(os.path.join(ROOT, 'results/summary.tsv'), 'w') as f:
    f.write('\t'.join(keys) + '\n')
    for r in rows:
        f.write('\t'.join('' if r.get(k) is None else ('%.3f' % r[k] if isinstance(r[k], float) else str(r[k])) for k in keys) + '\n')
print(open(os.path.join(ROOT, 'results/summary.tsv')).read())
