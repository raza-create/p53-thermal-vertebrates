import csv
def read_fa(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith('>'):
            n = l[1:].split('|')[1] if '|' in l else l[1:].split()[0]; s[n] = []
        else: s[n].append(l)
    return {k: ''.join(v) for k, v in s.items()}
aln = read_fa('vert_p53_aln.fasta')
meta = {}
with open('uniprot_vert_p53.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    lincol = [c for c in r.fieldnames if 'ineage' in c][0]
    for row in r:
        lin = row[lincol]
        cls = next((c for c in ['Mammalia','Aves','Amphibia','Actinopteri','Lepidosauria','Testudines','Crocodylia','Chondrichthyes'] if c in lin), 'Other')
        meta[row['Entry']] = (row['Organism'].split(' (')[0], cls)
H = aln['P04637']; pos2col = {}; p = 0
for c, a in enumerate(H):
    if a != '-': p += 1; pos2col[p] = c
dbd = list(range(94, 293)); key = [175,176,179,220,238,242,245,248,249,273,282]
rows = []; dbd_fa = open('vert_p53_DBD.fasta','w'); kp = open('key_positions.tsv','w')
kp.write('accession\tspecies\tclass\t' + '\t'.join(f'h{k}' for k in key) + '\n')
for acc, seq in aln.items():
    sp, cls = meta.get(acc, (acc, 'NA'))
    own = {}; q = 0
    for c, a in enumerate(seq):
        if a != '-': q += 1; own[c] = q
    cols = [pos2col[i] for i in dbd]
    ng = [c for c in cols if seq[c] != '-']
    cov = len(ng) / len(cols)
    idt = sum(seq[c] == H[c] for c in ng) / max(1, len(ng))
    keyres = ' '.join(f"{H[pos2col[k]]}{k}{seq[pos2col[k]]}" for k in key)
    rows.append((cls, sp, acc, len(seq.replace('-','')), cov, idt, keyres))
    dbd_fa.write(f">{acc}|{sp.replace(' ','_')}|{cls}\n{''.join(seq[c] for c in range(pos2col[94], pos2col[292]+1)).replace('-','')}\n")
    kp.write(f"{acc}\t{sp}\t{cls}\t" + '\t'.join(f"{seq[pos2col[k]]}{own.get(pos2col[k],'-')}" for k in key) + '\n')
dbd_fa.close(); kp.close()
order = {'Mammalia':0,'Aves':1,'Amphibia':2,'Actinopteri':3}
print(f"{'class':<12}{'species':<28}{'acc':<8}{'len':>4}{'DBDcov':>8}{'DBDid':>7}  key residues (human pos -> species aa)")
for r_ in sorted(rows, key=lambda x: (order.get(x[0],9), -x[5])):
    print(f"{r_[0]:<12}{r_[1][:27]:<28}{r_[2]:<8}{r_[3]:>4}{r_[4]:>8.2f}{r_[5]:>7.2f}  {r_[6]}")
