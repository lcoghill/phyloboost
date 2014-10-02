from Bio import SeqIO
import MySQLdb as mdb
import glob
import re



def get_ti(gi, db, gi_taxid_table) :

    sql = "".join(["SELECT taxid FROM ", gi_taxid_table, " WHERE gi=", gi])
    db.execute(sql)
    record = db.fetchone()
    
    if record:
        ti = str(record[0])
    else:
        ti = 0

    return ti

def get_genbank_data(gi, db, eukaroyte_table) :

    sql = "".join(["SELECT organism, acc, mol_type, division FROM ", eukaroyte_table, " WHERE gi=", gi])
    db.execute(sql)  # execute the above sql query
    record = db.fetchone()  # fetch all records that meet the query criteria and place them in the list record
    
    if record:
        tax_name = str(record[0])
        acc = str(record[1])
        mol_type = str(record[2])
        division = str(record[3])

    else:
        tax_name = "Unknown"
        acc = "Unknown"
        mol_type = "Unknown"
        division = "Unknown"
    
    return tax_name, acc, mol_type, division

def estimate_genome(gi, db, eukaroyte_table) :
    sql = "".join(["SELECT description FROM ", eukaroyte_table, " WHERE gi=", gi])
    db.execute(sql)
    record = db.fetchone()
    gene = get_gene(acc, db, genes_table)

    if record:
        description = str(record[0]).replace(",", "")
        mito = ["mito", "cyto"]
        nucl = ["nucl", "intron", "elongation", "exon", "EF-1"]
        chloro = ["chloro", "photosy"]
        if any(word in description for word in mito) :
            genome = "mitochondrial"

        elif any(word in description for word in nucl) :
            genome = "nuclear"

        elif any(word in description for word in chloro) :
            genome = "chloroplast"

        else:
            genome = "Unknown"

        if gene is "Unknown" :

            des_list = description.split(" ")
            for word in des_list:
                if "(" and ")" in word :
                    gene = word.replace("(", "").strip(")").lower()
                    if " " in gene :
                        gene = gene.split()[1]        

    else:
        descrption = "Unknown"
        genome = "Unknown"

    return description, genome, gene


def get_gene(acc, db, genes_table) :
    
    sql = "".join(["SELECT symbol, geneid FROM ", genes_table, " WHERE gi=", gi])
    db.execute(sql)
    record = db.fetchone()
    if record:
        gene = str(record[0]).lower()
    else:
        gene = "Unknown"
        
    return gene

def correct_taxid(ti, ti_dict) :

    if ti in ti_dict :
        new_ti = ti_dict[ti]
        return new_ti
    else:
        return ti



input_dir = ''
ti_error_file = 'merged.dmp'
out_file = 'phyloboost.csv'
db_host = 'localhost'
db_database = ''
db_user = ''
db_pass = ''
eukaroyte_table = ''
gi_taxid_table = ''
genes_table = ''
pb_version = "1.0"


## connect to the database and set a cursor
database = mdb.connect(host = db_host,  # your host, usually localhost
                     user = db_user,  # your username
                     passwd = db_pass,  # your password
                     db = db_database)  # name of the data base
db = database.cursor()  # define the db cursor


## get list of all cluster files in directory
cluster_files = glob.glob(input_dir + "*.fasta")
out_handle = open(out_file, 'w')

## builds a dict of ncbi deprecated ti values and their replacements
ti_handle = open(ti_error_file, 'r')
ti_dict = {}
for line in ti_handle :
    ti_dict[line.split("|")[0].replace("\t", "")] = line.split("|")[1].replace("\t", '')



for f in cluster_files :
    print "Converting file %s..." %f
    
    in_handle = open(f, 'r')
    
    for record in SeqIO.parse(in_handle, 'fasta') :
        ## for each sequence in cluster files
        ci = f.split("_ci")[1].replace(".fasta", "")
        gi = record.id.split("_")[0].replace("gi", "")
        ti = get_ti(gi, db, gi_taxid_table)
        ti = correct_taxid(ti, ti_dict)
        tax_name, acc, mol_type, division = get_genbank_data(gi, db, eukaroyte_table)
        sequence = str(record.seq)
        description, genome, gene = estimate_genome(gi, db, eukaroyte_table)
        out_string = ",".join([ti,tax_name,ci,sequence,genome,gene,acc,gi,pb_version,mol_type,division,description])
        out_handle.write(out_string+"\n")
    in_handle.close()

out_handle.close()
