import os
import re
import glob
import itertools
import multiprocessing
from multiprocessing import Process
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

###### FUNCTIONS ######

def parse_blastxml(file_name, results_path):
    # open blast output file	
    print "Parsing blast file %s" % file_name
    blast_list = []
    gi_list = []
    seq_list = []
    # iterate through the file
    for record in NCBIXML.parse(open(file_name)) :
        if record.alignments :
            for align in record.alignments :
                for hsp in align.hsps :
                    if align.length <= 7500: 
                        id_string = align.hit_id.split("gi|")
                        id_string = "".join(id_string)
                        id_string = id_string.encode("ascii")
                        id_string = id_string.split("|")
                        gi_list.append(id_string[0])

    return gi_list

## Perform the local blast searches
def blast_search(file_names, results_path, input_path):
	count = 0
	for f in file_names:
		print "Starting Blastn Search against %s....." %f
		fasta_file = "".join([input_path,"fasta/",f,".fasta"])
		gi_list = []
		## Blast using the FASTA set, and assing the XML string result to result_handle
		blast_file_id = "".join([results_path,"blast-results/", f, ".xml"])
		blastn_cline = NcbiblastnCommandline(query=fasta_file, db=blast_db, evalue=10, outfmt=5, out=blast_file_id)
		print "Writing results to a Genbank XML file....."
		## write the XML string to a file for later parsing
		stdout, stderr = blastn_cline()
		print "Blast search complete."
		## call function to parse xml files for list of gi values of interest
		print ("Getting list of gi values from blast results..."),
		gi_list = parse_blastxml(blast_file_id, results_path)
		## Write gi list file for each search
		gi_file = "".join([results_path,"blast-results/",f,".txt"])
    		with open(gi_file, 'w+') as gi_file:
        		for g in gi_list:
            			gi_file.write("%s\n" % g) # call to write the filtered fasta list to the correct fasta file
		gi_file.close()
		print "Complete"
		## Delete the blast xml file
		os.remove(blast_file_id)
		count = count + 1


def mygrouper(n, iterable):
	args = [iter(iterable)] * n
	return ([e for e in t if e != None] for t in itertools.izip_longest(*args))

def parse_config():
    config_params = [];
    config = open('private/config', 'r') 
    for line in (line for line in config if not line.startswith('###')):
        line = line.rstrip('\n')
        line = line.split("=")
        config_params.append (line[1])

    return config_params

def get_filenames(file_name):
    with open (file_name, "r") as fasta_file:
    ## parse input to get a clean file_id
        file_id = re.search('fasta/(.+?).fasta', file_name).group(1)
        return file_id

###


## Gather import configuration information from the config file
params = parse_config() # retreive the params from the config file

## Retreive NCBI-BLAST NT Database location from Config file
blast_db = params[4]


## Get a count of the number of fasta files in the appropriate directory

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("".join([params[7],"fasta/*.fasta"]))
file_count = float(len(fasta_files))
print "%s files successfully found.\n" %file_count
clean_filenames = []
for f in fasta_files:
    clean_filenames.append(get_filenames(f))

## Parallel Processing Code
core_count = multiprocessing.cpu_count()
print "There are %s cores on this system." %core_count
core_number = raw_input("Number of cores to use(1): ") # Ask the user how many cores to use

if not core_number:
	core_number = "1"

core_number = float(core_number)
print "Using %s cores." %core_number
files_per_core = int(round(file_count / core_number))

fasta_list = list(mygrouper(files_per_core, range(int(file_count)))) #create a matrix equal to the number of fasta_files per core with a total length of the number of fasta_files



## asign the file names to the proper elements inside of the fasta_list
track_fasta = 0
core_count = 0
for l in fasta_list:
	count = 0
	while count < len(fasta_list[core_count]):
		fasta_list[core_count][count] = clean_filenames[track_fasta]
		track_fasta = track_fasta + 1
		count = count + 1
	
	core_count = core_count + 1
## Split job into appropriate number of cores
processes = []
for c in xrange(len(fasta_list)):
	p = Process(target=blast_search, args=(fasta_list[c],params[6],params[7]))
	p.start()
	processes.append(p)

for p in processes:
	p.join()
