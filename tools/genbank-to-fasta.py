from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import sqlite3



### this script accepts a genbank(full) formatted file and gi_to_ti.db
### sqlite database containing the gi_ti dump from NCBI.
### it will parse the genbank file, reconstructing a fasta record properly
### formatted for the phyloboost pipeline.
### it is also an additional check to remove any non-eukaryotic DNA sequences.


genbank_file = 'genbank.gb' # sequences from Genbank in 'genbank (full)' format
fasta_file = 'eukaryotic_genbank.fas' # out file to save records in fasta format
sqlite_db = 'gi_to_ti.db'


con = sqlite3.connect(sqlite_db)
cur = con.cursor()
fasta_handle = open(fasta_file, 'a')
total_recs = 1
total_converted = 0
for record in SeqIO.parse(open(genbank_file), 'genbank') :
    total_recs += 1
    print record.annotations
    gi = record.annotations['gi']
    accession = record.id
    description = record.description
    sequence = record.seq
    cur.execute('SELECT * FROM gi_to_ti WHERE gi=%s' %int(gi))
    result = cur.fetchone()
    if result and 'RNA' not in description :
        ti = result[1]
        seq_id = "gi" + gi + "_ti" + str(ti) + "_acc" + accession
        seq_record = SeqRecord(Seq(str(sequence)), id=seq_id, description = description)
        SeqIO.write(seq_record, fasta_handle, 'fasta')
        print seq_record
        total_converted += 1

fasta_handle.close()
print "A total of %i / %i records were converted to FASTA format." %(total_converted, total_recs)
