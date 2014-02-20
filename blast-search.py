import multiprocessing
from multiprocessing import Process
import os
import re
import glob
import itertools
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import time



## Perform the local blast searches
def blast_search(file_names):
	count = 0
	for f in file_names:
		print "Reading original FASTA file %s....." %f[6:],
		print "Complete.\n"
		with open (f, "r") as fasta_file:
			## parse input to get a clean file_id
			sequences = fasta_file.read()
			fasta_file.close()
			print "Starting BLASTN Search against %s....." %f[6:],
			## Blast using the FASTA set, and assing the XML string result to result_handle
			blast_file_id = "".join(["blast-results/", f[6:-6], ".xml"])
			blastn_cline = NcbiblastnCommandline(query=f, db=blast_db, evalue=10, outfmt=5, out=blast_file_id)
			print "Complete."
			print "Writing results to a BLAST.XML file.....",
			## write the XML string to a file for later parsing
			stdout, stderr = blastn_cline()
			print "Complete."
			count = count + 1


def mygrouper(n, iterable):
	args = [iter(iterable)] * n
	return ([e for e in t if e != None] for t in itertools.izip_longest(*args))


## Retreive NCBI-BLAST NT Database location from Config file
config = [];
db_params = open('private/config', 'r')

for line in db_params:
    line = line.rstrip('\n')
    line = line.split("=")
    config.append (line[1])

blast_db = config[4]


## Get a count of the number of fasta files in the appropriate directory

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = float(len(fasta_files))
print "%s files successfully found.\n" %file_count



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
		fasta_list[core_count][count] = fasta_files[track_fasta]
		track_fasta = track_fasta + 1
		count = count + 1
	
	core_count = core_count + 1

## Split job into appropriate number of cores
processes = []
for c in xrange(len(fasta_list)):
	p = Process(target=blast_search, args=(fasta_list[c],))
	p.start()
	processes.append(p)

for p in processes:
	p.join()
