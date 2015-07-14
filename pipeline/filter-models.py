import sqlite3
from Bio import SeqIO



con = sqlite3.connect('gi_taxid_nucl.db')
def gi2ti(gi):
    q = 'select ti from gi_to_ti where gi = {}'.format(gi)
    try:
        return con.execute(q).fetchone()[0]
    except:
        return ''



## get models
merg_handle = open('merged.dmp', 'r')
merged = {}
print "Getting merged ti values..."
for line in merg_handle :
    rec = line.strip().split("\t|\t")
    merged[int(rec[0])] = int(rec[-1].replace("\t|", ""))

mhandle = open('model_organisms.txt', 'r')
models = []
print "Collecting and correcting model organism ti values..."
count = 0
for line in mhandle :
    ti = int(line.strip())
    if ti in merged :
        ti = merged[ti]
        count += 1
    models.append(int(line.strip()))

print "%i model organisms corrected." %count

print "Filtering model organisms from fasta..."
handle = open('eukaryotic_genbank.fas', 'r')

total = 0
kept = 0
for record in SeqIO.parse(handle, 'fasta') :
    ti = int(record.id.split("_")[-1][2:])
    total += 1
    out = open('expanded-fasta-all-no-models.fas', 'a')
    if ti not in models :
        SeqIO.write(record, out, 'fasta')
        kept += 1


print "A total of %i / %i records were copied." %(kept, total)
print "A total of $i records were discarded as being model organisms" %(total - kept)
out.close()
