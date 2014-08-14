from Bio import SeqIO
import glob
import MySQLdb
import re

def est_genome(description) :
    
    if "mito" in description:
        genome = "mitochondrial"

    elif "nucl" in description:
        genome = "nuclear"

    elif "chloro" in description:
        genome = "chloroplast"

    else:
        genome = "ignotus"

    return genome

def get_acc(gi, conn, cur) :

    cur.execute("SELECT acc FROM acc_to_gi WHERE gi=%s" %gi)
    row = cur.fetchone()

    if row:
        return row[0]

    else:
        return "ignotus"

def get_locus(gi, conn, cur) :

    cur.execute("SELECT geneid, symbol FROM genes WHERE gi=%s" %gi)
    row = cur.fetchone()

    if row:
        locus = row[1]
        geneid = row[0]
        return (locus, geneid)

    else:
        locus = "ignotus"
        geneid = "00000"

        return (locus, geneid)

def get_genome(gi, conn, cur) :

    cur.execute("SELECT description FROM eukaryotes WHERE gi=%s" %gi)
    row = cur.fetchone()

    if row:
        genome = est_genome(str(row[0]))
        return genome

    else:
        genome = "ignotus"
        return genome



conn = MySQLdb.connect(user="root", passwd="reelab14", db="genbank")
cur = conn.cursor()

fasta_dir = '/home/lcoghill/Dev/phyloboost/raw-fasta/'
fasta_files = glob.glob(fasta_dir + "*.fasta")
genome_f = open('genome_estimates.csv', 'a')

for fasta in fasta_files:

    handle = open(fasta, 'r')
    print "Estimating genomes for file %s..." %fasta
    for record in SeqIO.parse(handle, "fasta") :
        gi = int(record.id[3:])
        ci = int(re.search('_ci(.*).fasta', fasta).group(1))
        acc = get_acc(gi, conn, cur)
        locus, geneid = get_locus(gi, conn, cur)
        genome = get_genome(gi, conn, cur)
        row = ",".join([str(gi), str(ci), acc, locus, str(geneid), genome])
        genome_f.write(row + '\n')
    
    handle.close()

genome_f.close()

