import os
import re
import sets
import glob
import logging
import MySQLdb
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline

## Retreive NCBI-BLAST NT Database location from Config file
config = [];
db_params = open('private/config', 'r')

for line in db_params:
    line = line.rstrip('\n')
    line = line.split("=", 1)
    config.append (line[1])

blast_db = config[4]



print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successfully found.\n" %file_count

## WARN IF THERE ARE A LOT OF FASTA FILES
if file_count > 2:
	print "WARNING:\n"
	print "There are %s FASTA files to be processed." %file_count
	print "This process can take quite a long time...\n\n"

### Start processing BLAST searches if the user agrees
count = 0
for f in fasta_files:
	print "Reading original FASTA file %s" %f[6:]
	print "Successfully imported.\n"
	with open (f, "r") as fasta_file:
		## parse input to get a clean file_id
		print "Step %s / %s beginning." % (count + 1, file_count)
		sequences = fasta_file.read()
		fasta_file.close()
		print "Starting BLASTN Search against %s" %f[6:]
		## Blast using the FASTA set, and assing the XML string result to result_handle
		blast_file_id = "".join(["blast-results/", f[6:-6], ".xml"])
		blastn_cline = NcbiblastnCommandline(query=f, db=blast_db, evalue=10, outfmt=5, out=blast_file_id)
		print "Done."
		print "Writing results to a BLAST.XML file..."
		## write the XML string to a file for later parsing
		stdout, stderr = blastn_cline()
		print "Step %s / %s complete." % (count + 1, file_count)
		count = count + 1


	print "BLAST search successfully completed.\n"

print "Process Complete."

