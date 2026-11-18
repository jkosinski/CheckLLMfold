# Fable 5.2 contact prediction on proteins with and without PDB relatives

Date: 2026-11-18. Protocol: Sergey Ovchinnikov's prompt (shared chat claude.ai/share/068377d7-bace-48bf-8102-c0f144da3887),
text in `../prompt_template.md`. X thread: https://x.com/jankosinski/status/2122997527192416440

## Setup

- Model: `claude-fable-5-2`, effort high, run headless with `claude -p` (`../scripts/run_fable.sh`).
- Tools: Bash/Read/Write only, OS sandbox with no network, no MCP servers, user CLAUDE.md excluded.
- Working directory `runs/jobs/<random id>/work`, so the path does not reveal the protein. The first ubiquitin run used a path
  containing "ubiquitin_1UBQ"; it was discarded and is not in this repository.
- Changes to Sergey's prompt: the py2Dmol artifact step is dropped; the final JSON also returns the final contact list
  ("final_contacts") so the model can be rebuilt with the same `fold()` code. The number of ranked pairs is 286 for
  ubiquitin, as in his prompt, and round(286/76 x L) for the others.
- Scoring (`../scripts/evaluate.py`): contact = C-beta (C-alpha for Gly) < 8 A. Precision of the ranked pairs at top L/5,
  L/2, L for |i-j| >= 6 and for |i-j| >= 24. TM-score of the CA model rebuilt from the final contacts, for both mirror
  images, because `fix_hand()` in the protocol picks the wrong hand for some proteins (cyt b562, thioredoxin, and tenascin
  even with the true contacts). Ceiling = `fold()` on the true contacts.

## Targets

| Target | L | Class | Reference | In PDB |
|---|---|---|---|---|
| Ubiquitin | 76 | α/β | 1UBQ chain A (crystal) | yes, famous |
| Cytochrome b562 | 106 | all-α | 256B chain A (crystal) | yes, famous |
| Tenascin FN3 domain | 89 | all-β | 1TEN chain A (crystal) | yes, famous |
| Thioredoxin (E. coli) | 108 | α/β | 2TRX chain A (crystal) | yes, famous |
| A0A916UKZ0, DUF2489 (PF10675), *Parundibacterium terreum* | 108 | all-α | AFDB v6, pLDDT 97.4, mean PAE 2.4 A | no |
| A0A8J2YCU1, DUF2958 (PF11171) | 89 | all-β | AFDB v6, pLDDT 91.1, mean PAE 3.8 A | no |
| A0A5D0P002, DUF4288 (PF14119), *Actinomadura chibensis* | 107 | α/β | AFDB v6, pLDDT 94.4, mean PAE 3.6 A | no |

How "not in PDB" was established for the three AFDB proteins (`../data/candidates/`):
- RCSB MMseqs2 sequence search against the PDB: no hit at E < 1.
- Pfam family has 0 structures in InterPro and belongs to no clan.
- Foldseek (3Di+AA) against PDB100: no hit with E < 0.5 (best E 1.1, 0.52, 0.59).
- Foldseek TM-align mode: best TM-score normalised by query length 0.60, 0.52, 0.54. These are shape-level matches to
  unrelated proteins (other helical bundles, an Sm-like barrel, a ribosomal protein), not homologs.

Candidates came from 500 small bacterial DUF proteins in UniProt with an AFDB model and no PDB cross-reference, filtered to
pLDDT >= 90. Of 27 checked with Foldseek, most had a structural homolog in the PDB; DUF559 was dropped because its Pfam
family has 2 structures. The Foldseek server rate-limited the remaining queries.

## Results

`summary.tsv` has all metrics; `fig_summary.png` and `fig_contact_maps.png` show them.

| Target | Runs | Top-L, \|i-j\| >= 6 | Top-L, \|i-j\| >= 24 | TM, better hand (ceiling) |
|---|---|---|---|---|
| Ubiquitin | 3 | 0.62, 0.78, 0.62 | 0.45, 0.68, 0.51 | 0.44, 0.52, 0.44 (0.79) |
| Ubiquitin, Sergey's chat | 1 | 0.80 | 0.66 | n/a |
| Cytochrome b562 | 1 | 0.24 | 0.18 | 0.42 (0.81) |
| DUF2489, all-α | 3 | 0.06, 0.06, 0.04 | 0.04, 0.05, 0.06 | 0.21, 0.22, 0.22 (0.79) |
| Tenascin FN3 | 1 | 0.45 | 0.47 | 0.33 (0.79) |
| DUF2958, all-β | 2 | 0.37, 0.43 | 0.04, 0.07 | 0.20, 0.20 (0.73) |
| Thioredoxin | 1 | 0.52 | 0.49 | 0.47 (0.84) |
| DUF4288, α/β | 3 | 0.24, 0.30, 0.44 | 0.20, 0.07, 0.33 | 0.24, 0.29, 0.30 (0.61) |

Chance precision at |i-j| >= 6 is 0.02-0.06. The "Sergey's chat" row scores the 286 pairs from his shared chat; that
output is not included here. To rescore it, save his final answer as `runs/sergey_shared_chat/final.txt` and run
`scripts/aggregate.py`.

- The control reproduces: on ubiquitin Fable gets the secondary structure and the strand registers right, with top-L/5
  precision 0.80-1.00.
- On the three proteins with no PDB relative, long-range precision is at chance for the all-α and all-β proteins and
  0.07-0.33 for the α/β one. No model reaches the fold (TM 0.20-0.30).
- Size- and class-matched famous proteins do much better at long range: 0.47 (tenascin) against 0.04-0.07 (DUF2958), 0.49
  (thioredoxin) against 0.07-0.33 (DUF4288), 0.18 (cyt b562) against 0.04-0.06 (DUF2489).
- Contacts that follow from local sequence patterns survive on the unknown proteins: helical i/i+3, i/i+4 contacts on the
  all-α protein (precision 0.78-0.87), and β-hairpins joined by short turns. Secondary structure is sometimes misread
  (a helix called two strands in DUF2489; helices predicted in the all-β DUF2958).

## What Fable knows about the sequences

### Strand labels in the contact runs

Fable's visible text in the contact runs never names a protein, a family or a PDB entry. Its code (tool calls in
`runs/jobs/*/transcript.jsonl`) does show fold knowledge on two of the famous proteins:

- Tenascin (`job_3080fa44`): the seven strands are labelled A, B, C, C', E, F, G, with no D. These are the standard FN3
  strand names. The pairings it sets, A-B, B-E, C-C', C-F and F-G, are the two FN3 sheets: all five are at 4.1-4.3 A
  minimum CA-CA distance in 1TEN, and the other strand combinations at 6.2 A or more. The write-up says the strands were
  read "from the alternation of hydrophobic residues".
- Thioredoxin (`job_47423b69`): b1-b3 and b2-b3 parallel, b2-b4 antiparallel, b4-b5 hairpin. This is the thioredoxin sheet
  (4.3-4.6 A in 2TRX; the alternatives 8.9 A or more). The write-up presents the secondary structure as read from the
  sequence.
- The three AFDB proteins get plain sequential labels (b1, b2, a1, ...; A-E with D included).

The thinking blocks cannot be used for this: in 7 of the 14 runs their text is empty, in the other 7 it is a one- or
two-sentence summary of a round.

### Recognition probe

Fable 5.2 got each full sequence, with no tools, and was asked for the protein name, accession, organism, Pfam family,
fold and PDB IDs (`../scripts/recognize.py`, same headless setup as `nextaa/scripts/probe.py`). Questions are in
`../data/recognition/`, replies and `summary.tsv` in `recognition/`, scoring in `../scripts/recognition_score.py`.

Contact-test proteins, 3 questions at effort low and 3 at effort high:

| Protein | Accession | Pfam | Organism | PDB |
|---|---|---|---|---|
| Ubiquitin, cyt b562, tenascin FN3, thioredoxin | 6/6 each | 6/6 each | 6/6 each | 6/6 each |
| A0A916UKZ0 (DUF2489) | 0/6 | 2/6 | 0/6 | - |
| A0A8J2YCU1 (DUF2958) | 0/6 | 6/6 | 0/6 | - |
| A0A5D0P002 (DUF4288) | 0/6 | 6/6 | 0/6 | - |

For the famous proteins the PDB IDs include the reference entry used here (1UBQ, 256B, 1TEN, 2TRX) in all 24 answers. For
the DUF proteins Fable says the family is assigned "from recalled family motifs, not a database search".

Proteins from the next-residue test, effort low, one question each (`nextaa/README.md`):

| Set | Proteins | Accession right | Pfam family right (chance) |
|---|---|---|---|
| Well-known (ubiquitin, GFP, p53, ...) | 12 | 12 (protein named in all) | - |
| Pre-cutoff, UniRef50 cluster > 50 | 40 | 0 | 75% (3%) |
| Pre-cutoff, no UniRef50 relative | 40 | 0 | 62% (2%) |
| Post-cutoff (Tu & Chen's D set) | 56 | 0 | 75% (4%) |

Chance: the same score with the answers shuffled among the proteins of a set. Pfam counts as right if the true family
(accession or name, e.g. DUF222) is among the families given. Fable gave no accession for any of these proteins and
never the right species; the genus was right only in guesses, mostly from composition ("a high-GC actinobacterium such
as Streptomyces"), in 18%, 10% and 2% of the three sets.

Reading:
- Fable identifies the famous proteins down to the UniProt and PDB entry, so their contact maps can draw on known
  structures.
- For obscure proteins it names the Pfam family about as often for proteins first public after its training cutoff as for
  older ones, and never gives an accession. This looks like family-level knowledge, similar to a profile search, and not
  recall of individual sequences.
- It named the family of two of the three AFDB proteins every time, yet contact precision on them stayed near chance
  (DUF2958) or low (DUF4288). These families have no structure in the PDB. Naming the family did not give good contacts
  without a known structure of the family.
- In the next-residue test, the first-residue gap (real minus shuffled) does not follow family recognition in a consistent
  way: pre-cutoff 0.28 when the family was named and 0.20 when not, post-cutoff 0.08 and 0.26 (14-42 proteins per group).

## Confidence of the AFDB references

| Target | Mean pLDDT | Residues >= 90 | Residues < 70 | Mean PAE | PAE over native contacts \|i-j\| >= 24 |
|---|---|---|---|---|---|
| DUF2489 | 97.5 | 99% | none | 2.4 A | 1.1 A |
| DUF2958 | 91.1 | 72% | 1, 88, 89 | 3.8 A | 1.8 A |
| DUF4288 | 94.4 | 92% | 1-2, 106-107 | 3.6 A | 1.3 A |

Re-scoring with only pairs whose residues both have pLDDT >= 90 and PAE < 5 A in both directions gives long-range top-L
precision 0.04-0.06 (DUF2489), 0.07-0.08 (DUF2958, chance 0.066 on this subset) and 0.08-0.33 (DUF4288), the same as with
all pairs.

## Figures for X

`x_thread/`, 1600x900, in thread order: x1 + x1b (ubiquitin contact maps and structures), x2 (selection), x3a/b/c
(class-matched contact maps) each with its x3*_structures image, x4 (summary). Built by `scripts/figures_x.py`; structure
renders by `scripts/prep_structures.py`, `scripts/model_ss.py` and ChimeraX (`x_thread/structures/render.cxc`).

Structure images: the reference is a normal cartoon from its full-atom coordinates. The model is the CA trace from `fold()`
on Fable's final contacts (better mirror image), superposed on the reference by TM-score. Its helices and strands are read
from Fable's own contacts (strand: ladder of contacts with |i-j| >= 6, split at hairpin turns; helix: three consecutive
i/i+3 and i/i+4 turns, >= 6 residues), not from the reference. For proteins with several runs the run with the best
long-range precision is shown, as in the contact maps.

## Caveats

- One replicate per famous control; three (two for DUF2958) per unknown protein.
- The unknown proteins are references from AlphaFold, not experiments, though with very high confidence.
- A famous protein may be easier for reasons other than familiarity. The class-matched pairs reduce this but do not remove
  it.
- Fable's reasoning is not visible: thinking blocks come back empty or as short summaries. See "What Fable knows about the
  sequences" for what its code and a direct question show.
