from Bio import SeqIO
import sqlite3
import glob



expanded_fasta_dir = ''
db = sqlite3.connect('gi_taxid_nucl.db')
cursor = db.cursor()
## get list of all fasta files
fasta_files = glob.glob(expanded_fasta_dir+"*.fas")



count = 1
for fasta in fasta_files :
    print "Renaming all records in file %s / %s..." %(count, len(fasta_files))
    new_records = []    
    ## open file in read format
    records = list(SeqIO.parse(fasta, "fasta"))
    
    ## get list of all records
    for rec in records :
        
        ## find id in sqlite to get TI
        cursor.execute('''SELECT ti FROM gi_to_ti WHERE gi=%s''' %rec.id[3:])
        row = cursor.fetchone()
        if row :
            ti = row[0]
            rec.id = "gi"+str(rec.id[3:])+"_ti"+str(ti)
            rec.name = ""
            rec.description = ""
            print rec.id
            new_records.append(rec)

    handle = open(fasta, "w")
    SeqIO.write(new_records, handle, "fasta")
    handle.close()
    count += 1

print "Rename complete."
