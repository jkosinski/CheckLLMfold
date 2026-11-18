"""UniProtKB members of a Pfam family with creation date. Usage: uniprot_family.py PFAM NAME OUTDIR
-> OUTDIR/NAME.fasta, OUTDIR/NAME.dates.tsv (accession, createdate as YYYY/MM/DD, length, title)"""
import sys, os, urllib.request, urllib.parse
pf, name, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
url = 'https://rest.uniprot.org/uniprotkb/stream?' + urllib.parse.urlencode(
    {'query': f'xref:pfam-{pf}', 'format': 'tsv', 'fields': 'accession,date_created,length,protein_name,sequence'})
txt = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'python-urllib/3 research-script'}), timeout=600).read().decode()
lines = txt.strip().split('\n')[1:]
with open(os.path.join(out, name + '.fasta'), 'w') as fa, open(os.path.join(out, name + '.dates.tsv'), 'w') as dt:
    dt.write('accession\tcreatedate\tlength\ttitle\n')
    for l in lines:
        acc, dc, ln, pn, seq = l.split('\t')
        fa.write(f'>{acc}\n{seq}\n'); dt.write(f'{acc}\t{dc.replace("-", "/")}\t{ln}\t{pn}\n')
print(name, pf, len(lines))
