import os
import re
import glob
import MySQLdb as mdb
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Entrez



def get_filenames(file_name):
    ## parse input to get a clean file_id
    file_id = "ti" + re.search('/ti(.+?).txt', file_name).group(1)
    return file_id

def get_done_files(file_name):
    ## parse input to get a clean file_id
    file_id = "ti" + re.search('/ti(.+?).fasta', file_name).group(1)
    return file_id

def fetchseq(gis, entrez_email):
    ## function to fetch missing sequences and ids from web via entrez
    ## not currently in use
    missing_seqs = []
    missing_tis = []
    Entrez.email = entrez_email
    ## Upload the list of IDs 
    try:
        request = Entrez.epost("nucleotide",id=",".join(gis))
        result = Entrez.read(request)
        webEnv = result["WebEnv"]
        queryKey = result["QueryKey"]
        handle = Entrez.efetch(db="nucleotide", retmode="xml", webenv=webEnv, query_key=queryKey)
        for r in Entrez.parse(handle):
            missing_seqs.append(r["GBSeq_sequence"])
    except:
        pass
    return missing_seqs



## Gather import configuration information from the config file
db_host = "localhost" # database host usually 'localhost'
db_user = "" # database user
db_passwd = "" # database password for db_user
db_name = "" # database to use

fasta_dir = "" # fasta files input directory
gi_dir = "" # gi files input directory
done_dir = "" # check this directory to skip files already done
out_dir = "" # location to store output files
entrez_email = "" # email address used for entrez web searches


## Get List of FASTA FIles
fasta_files = glob.glob("".join([fasta_dir,"*.fasta"]))
## Get List of GI FIles
gi_files = glob.glob("".join([gi_dir,"*.txt"]))
## Get List of Already Expanded Files
done_files = glob.glob("".join([done_dir,"*.fasta"]))
## Get List of Clean ID Labels
clean_filenames = []
for f in gi_files:
    clean_filenames.append(get_filenames(f))

done_filenames = []
for f in done_files:
    done_filenames.append(get_done_files(f))

## Open connection to the database
database = mdb.connect(host=db_host,
                      user=db_user,
                      passwd=db_passwd,
                      db=db_name)

genbank_db = database.cursor() # define the db cursor
rec_count = 1

for c in clean_filenames:
    record_list = []
    if c not in done_filenames:
        print "Processing file %s, number %s / %s" %(c, rec_count, (len(clean_filenames)-len(done_filenames)))
        missing_gis = []
        missing_seqs = []
        ## Pull in original fasta to list 
        old_fasta = []
        new_fasta = []
        old_file = "".join([fasta_dir,c,".fasta"])
        for record in SeqIO.parse(old_file, "fasta"):
            gi = record.id.split("|")[1].strip('\n')
            fetch_ti_sql = "".join(["SELECT taxid FROM gi_to_taxid WHERE gi=",gi,";"])
            genbank_db.execute(fetch_ti_sql) #execute the above sql query
            ti = genbank_db.fetchone()
            if ti:
                ti = str(ti[0])
                sequence = str(record.seq)
                id = "".join(["gi",gi,"_ti",ti])
                record_list.append(SeqRecord(Seq(sequence[1:], IUPAC.unambiguous_dna), id=id, description=""))

        gi_list = []
        gi_file = "".join([gi_dir,c,".txt"])
        with open(gi_file) as f:
            gi_list = f.readlines()
            print "Retreiving sequences from local MySQL database..."
            for g in gi_list:
                ## Use GIs from List to Pull Seqs from MySQL
                fetch_seq_sql = "".join(["SELECT sequence FROM eukaryotes WHERE gi=",g.strip("\n"),";"])
                genbank_db.execute(fetch_seq_sql) #execute the above sql query
                sequence = genbank_db.fetchone()
                if sequence:
                    fetch_ti_sql = "".join(["SELECT taxid FROM gi_to_taxid WHERE gi=",g.strip("\n"),";"])
                    genbank_db.execute(fetch_ti_sql) #execute the above sql query
                    ti = genbank_db.fetchone()
                    if ti:
                        ti = str(ti[0])
                        sequence = str(sequence[0])
                        gi = "".join(["gi",g,"_ti",ti])
                        record_list.append(SeqRecord(Seq(sequence[1:], IUPAC.unambiguous_dna), id=gi, description=""))
                    else:
                        gi_list.remove(g)
                else:
                    missing_gis.append(g.strip("\n"))
        for r in record_list:
            r.seq = r.seq.lower()
        print "Writing new %s records to fasta file for cluster %s." %(len(record_list), c)
        fasta_file = "".join([out_dir,c,".fasta"])
        SeqIO.write(record_list, fasta_file, "fasta") #after done with all iterations, write the good seq record list to the same file
        rec_count = rec_count + 1
        print "-"*50
    else:
        print "File %s, already merged. Skipping." %c
        print "-"*50    
database.close()

print "Fasta-merge complete."
