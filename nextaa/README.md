# Next-amino-acid prediction by Fable 5.2

> **AI-generated.** Code, analysis and text by Claude Opus 5.5 in Claude Code; predictions by Claude Fable 5.2. Not
> thoroughly revised by a human. There may be mistakes.

Date: 2026-11-18. Xinming Tu and Leo Chen ([LLM-bioseq-probe](https://github.com/XinmingTu/LLM-bioseq-probe)) report
that Claude Fable 5.2 predicts the next amino acid of obscure bacterial proteins (Pfam DUF families) far above a
shuffled control. The question here: is this recall of sequences in public databases? X thread: https://x.com/jankosinski/status/2123006759358914813

## Setup

- Prompt: Tu & Chen's chat prompt, verbatim (`scripts/probe.py`). The model gets a prefix of the sequence and replies with
  the next 10 residues. Only the first residue is scored unless stated otherwise.
- Control: the same protein with its residues shuffled (seeded by accession), cut at the same position. Gap = real minus
  shuffled first-residue accuracy, averaged over Pfam families, with a 95% bootstrap CI over families (their statistic).
- Our runs: `claude-fable-5-2`, effort low, through headless Claude Code (`claude -p`) with no tools, no MCP servers,
  CLAUDE.md files excluded and an empty working directory. 1,428 calls, all replies in `results/*.jsonl`.
- Harness effect: headless Claude Code adds about 630 tokens of its own context and scores about 10 points below Tu &
  Chen's API runs on the same questions (`results/matched_summary.tsv`). Old and new proteins are only compared within
  one setup.
- Unseen proteins: Fable 5.2's training data ends in June 2026. Proteins first public after that come from RefSeq, which
  adds many new sequences every month. Some new RefSeq records carry old sequences, so each exact sequence was checked in
  UniParc. UniParc lags RefSeq by months.
- Relatives: number of UniProt proteins at 50% identity or more = UniRef50 cluster size minus one.

## Question sets

| Set | Proteins | File |
|---|---|---|
| Well-known proteins (ubiquitin, GFP, p53, KRAS, ...), 5 cuts each | 12 | `data/q_famous.tsv` |
| Pre-cutoff proteins with no UniRef50 relative, from Tu & Chen's main set | 116 | `data/q_pre_single.tsv` |
| Post-cutoff proteins (Tu & Chen's D set) | 56 | `data/q_post_D.tsv` |
| Divergence test: new proteins 70-95% identical to an older UniProtKB protein | 48 | `data/div_uniprot/q_div.tsv` |
| Same sites on the closest old relative | 48 relatives | `data/div_uniprot/q_rel.tsv` |

`data/div_run1/` is a first version of the divergence test with RefSeq relatives. Most of those were not in UniProtKB,
so it was replaced by the UniProtKB version (`results/div_fable52.jsonl` and `relself_fable52.jsonl` belong to it).

## Results

- Well-known proteins: right next residue at all 60 cuts, exact next 10 residues at 58. Shuffled: 10%.
- Tu & Chen's published Fable 5.2 replies, split by relatives in UniProt (`results/tu_reanalysis_by_uniref50.tsv`):

  | UniRef50 cluster size | Fable 5.2 gap | Opus 5.5 | ProGen3-3B | ProtGPT3-10B |
  |---|---|---|---|---|
  | > 50 | 54.1 | 36.0 | 47.1 | 42.6 |
  | 6-50 | 46.5 | 26.3 | 33.3 | 27.9 |
  | 2-5 | 36.1 | 20.5 | 18.2 | 12.7 |
  | 1 (no relative) | 36.4 | 22.4 | 23.8 | 18.9 |
  | post-cutoff | 22.1 | 13.4 | 12.0 | 6.7 |

  Post-cutoff gap as a fraction of the no-relative gap: Fable 0.61, Opus 0.60, ProGen3 0.51, ProtGPT3 0.35. Our own Fable
  runs on the same two sets: 24.6 and 15.0.
- Divergence test (`scripts/analyze_div_paired.py`, `scripts/analyze_div_neighbors.py`):
  - At sites where the new protein differs from its closest old relative, and Fable gives the right residue on the old
    relative itself (29 sites), Fable gives the old residue on the new protein 59% of the time and the new residue 7%.
  - At sites where the two agree and Fable is right on the old relative (23 sites), it is right on the new protein 96%.
  - The old residue it gives is also the residue of the 10 next-closest old relatives in 58% of cases, so it follows the
    subfamily, not one sequence.
  - First-residue accuracy is about the same on the new proteins (26%) and on their old relatives (27%).

Reading: Fable 5.2 behaves like a model trained on UniProt, best near known sequences and worse far from them, with
higher scores than the open protein language models tested. This does not separate memorization from what protein
language models do.

## Caveats

- Fable was run through headless Claude Code, not the API (no API key on the machine used).
- The divergence test is small (48 proteins, 192 sites), and no protein language model was run on it.
- New-protein status depends on RefSeq and UniParc dates at the time of the test (2026-11).

## Files

- `scripts/`: `probe.py` (runner), `make_questions.py`, `divergence_set.py`, `relative_questions.py`, analysis scripts
  (`stratify_tu.py`, `analyze_matched.py`, `analyze_div*.py`), data download helpers (`ncbi_family.py`,
  `uniprot_family.py`, `uniref_map.py`, `uniparc_first_seen.py`, `uniprot_status.py`), `figures_x.py`.
- `data/tu_outputs/`: replies published by Tu & Chen, copied from their repository (MIT licence, `LICENSE` there).
- `data/div_*/hits.m8`: MMseqs2 hits of new against old proteins, reduced to the rows the scripts read (hits of the
  test proteins with query coverage >= 0.5, top 50 per protein, plus each protein's closest relative). The RefSeq and
  UniProtKB family downloads (55 MB) are not included; `ncbi_family.py` and `uniprot_family.py` fetch them again, but
  the databases change over time. `old.fasta` and `new.fasta` hold the sequences that were searched.
- `results/x_thread/`: figures and text written for an X thread.

Run from this folder, e.g. `python3 scripts/probe.py data/q_famous.tsv results/famous_fable52.jsonl --jobs 6`. It uses
`../runs/sandbox_settings.json`; edit the CLAUDE.md path there first.
