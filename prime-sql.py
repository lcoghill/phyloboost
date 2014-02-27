import os
import MySQLdb as mdb
from Bio import SeqIO
from Bio import Entrez
import time
import mmap
import random
from collections import defaultdict



###### FUNCTIONS  #######
        
def parse_config():
    config_params = [];
    config = open('private/config', 'r') 
    for line in (line for line in config if not line.startswith('###')):
        line = line.rstrip('\n')
        line = line.split("=")
        config_params.append (line[1])

    return config_params

#####



## Gather import configuration information from the config file
params = parse_config() # retreive the params from the config file
print "Please enter the name of your nucleotide file from Genbank."
print "The script will look in %s (specified input location) for the file." % params[7]
file_name = raw_input("Filename: ")
file_name = "".join([params[7],file_name])
print "Opening %s" % file_name
print "-"*50
print "\n"


## Open connection to the database
database = mdb.connect(host=params[1],
                      user=params[2],
                      passwd=params[3],
                      db=params[0])

cluster_db = database.cursor() # define the db cursor

## Check of the table exists. If it does, drop it and re-create. Otherwise, just create the table.
cluster_db.execute("DROP TABLE IF EXISTS genbank_nucleotide")
cluster_db.execute("CREATE TABLE genbank_nucleotide(id INT PRIMARY KEY AUTO_INCREMENT, gi INT(20), accession VARCHAR(25), length int(10), title VARCHAR(25), sequence LONGTEXT)")

## Iterate through each record in the fasta file, parse out the important information and write that information to records in the mysql database.
f = open(file_name,'r+')
count = 0
for record in SeqIO.parse(file_name, "fasta") :
     seq_list = record.description.split("|")
     print seq_list
     cluster_db.execute('''INSERT into genbank_nucleotide (gi, accession, length, title, sequence) VALUES (%s, %s, %s, %s, %s)''', (seq_list[1], seq_list[3], len(record.seq), seq_list[4][1:], record.seq))
     database.commit()
     count = count + 1

f.close()

##disconnect from db
database.close()

print "-"*50
print "%s records successfully loaded into database." % count