import os
import re
import glob
import MySQLdb as mdb

def get_filenames(file_name):
    with open (file_name, "r") as fasta_file:
    ## parse input to get a clean file_id
        file_id = re.search('fasta/(.+?).fasta', file_name).group(1)
        return file_id

def parse_config():
    config_params = [];
    config = open('private/config', 'r') 
    for line in (line for line in config if not line.startswith('###')):
        line = line.rstrip('\n')
        line = line.split("=")
        config_params.append (line[1])

    return config_params


## Gather import configuration information from the config file
params = parse_config() # retreive the params from the config file
## Get List of FASTA FIles
fasta_files = glob.glob("".join([params[7],"fasta/*.fasta"]))
## Get List of GI FIles
gi_files = glob.glob("".join([params[6],"blast-results/*.txt"]))


## Get List of Clean ID Labels
clean_filenames = []
for f in fasta_files:
    clean_filenames.append(get_filenames(f))

## Open connection to the database
database = mdb.connect(host=params[1],
                      user=params[2],
                      passwd=params[3],
                      db=params[0])

db = database.cursor() # define the db cursor

for c in clean_filenames:
    ## Pull in original fasta to list 
    old_fasta = []
    new_fasta = []
    old_file = "".join([params[7],"/fasta/",c,".fasta"])
    with open(old_file) as f:
        old_fasta = f.readlines()
    ## Pull in list of GI values
    gi_list = []
    gi_file = "".join([params[6],"/blast-results/",c,".txt"])
    with open(gi_file) as f:
        gi_list = f.readlines()
        for g in gi_list:
            ## Use GIs from List to Pull Seqs from MySQL
            fetch_seq_sql = "".join(["SELECT sequence FROM genbank_nucleotide WHERE gi=",g,";"])
            new_fasta.append(g)
            new_fasta.append(db.execute(fetch_seq_sql)) #execute the above sql query
    print new_fasta 

### create new fasta list
### write to new_fasta file