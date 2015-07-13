from Bio import SeqIO
import glob


handle = open('outlier-results.txt', 'r')

outliers = {}
out_dir = ''

## collect all outlier results from find-outliers.py
for line in handle :
    rec = line.strip().replace("]", "").split(" [")
    tis = rec[1].replace(" ", "").split(",")
    outliers[rec[0].split("/")[-1]] = tis

## open the flagged alignment, remove any sequences that have the
## suspected misidentified TI values
for key, val in outliers.items() :
    sequences = list(SeqIO.parse(open(key), 'fasta'))
    good_seqs = []
    for s in sequences :
        ti = s.id.split("|")[1][2:]
        if ti not in val :
            good_seqs.append(s)

    name = out_dir + key
    SeqIO.write(good_seqs, out, 'fasta')
