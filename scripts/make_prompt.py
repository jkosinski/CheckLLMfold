"""Fill prompt_template.md for one sequence. Usage: make_prompt.py SEQUENCE NPAIRS OUT.md
Deviations from Sergey's prompt: the py2Dmol artifact line is dropped (no artifact tool in a headless run) and the
final JSON also carries the final-round contact list so the model can be rebuilt for scoring."""
import sys, re, pathlib

root = pathlib.Path(__file__).resolve().parent.parent
seq, npairs, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
t = (root / 'prompt_template.md').read_text()
t = t.replace('{SEQUENCE}', seq).replace('{NPAIRS}', str(npairs))
t = t.replace('{EXTRA}\n', 'Also give your final-round contact list (all separations, |i - j| >= 3) in the JSON as "final_contacts".\n\n')
t = t.replace('{EXTRA_JSON}', ', "final_contacts": [[<int>, <int>], ...]')
pathlib.Path(out).write_text(t)


def code_block():
    s = (root / 'prompt_template.md').read_text()
    return re.search(r'```python\n(.*?)```', s, re.S).group(1)
