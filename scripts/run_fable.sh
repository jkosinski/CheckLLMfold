#!/bin/bash
# Run Sergey's contact-prediction protocol headless with Claude Fable 5.2.
# Usage: run_fable.sh TARGET REP SEQUENCE NPAIRS [EFFORT]
# The model's working directory is runs/jobs/<random id>/work so the path does not reveal the target;
# runs/jobs.tsv maps job ids to targets.
# Sandbox: Bash/Read/Write only, no network, no MCP, user CLAUDE.md excluded, fresh empty working directory.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
T=$1; R=$2; SEQ=$3; N=$4; EFFORT=${5:-high}
J=job_$(openssl rand -hex 4)
D=$ROOT/runs/jobs/$J
mkdir -p "$D/work"
printf '%s\t%s\t%s\t%s\n' "$J" "$T" "$R" "$EFFORT" >> "$ROOT/runs/jobs.tsv"
python3 -I "$ROOT/scripts/make_prompt.py" "$SEQ" "$N" "$D/prompt.md"
cd "$D/work"
date -u +%FT%TZ > "$D/started"
claude -p --model claude-fable-5-2 --effort "$EFFORT" \
  --tools "Bash,Read,Write" --setting-sources project --strict-mcp-config \
  --settings "$ROOT/runs/sandbox_settings.json" \
  --system-prompt "You have a code tool: the Bash tool runs shell commands, with python3 and numpy installed. Use it to run code." \
  --permission-mode acceptEdits --no-session-persistence \
  --output-format stream-json --verbose < "$D/prompt.md" > "$D/transcript.jsonl" 2> "$D/stderr.txt"
date -u +%FT%TZ > "$D/finished"
python3 -I -c "
import json
for l in open('$D/transcript.jsonl'):
    try: o=json.loads(l)
    except Exception: continue
    if o.get('type')=='result':
        open('$D/final.txt','w').write(o.get('result') or '')
        print('$J $T rep$R', 'cost_usd', o.get('total_cost_usd'), 'turns', o.get('num_turns'), 'duration_s', (o.get('duration_ms') or 0)/1000)
"
