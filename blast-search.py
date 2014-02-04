import os
import re
import sets
import glob
import logging
import MySQLdb
import library.genbank
logging.basicConfig(level=logging.INFO)
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML

### A few import static variables
library.genbank.email = 'me@my.address.com'    

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successful found.\n" %file_count

print "WARNING:\n"
print "There are %s FASTA files to be processed." %file_count
print "This process can take quite a long time.\n\n"

for f in fasta_files:
	print "Reading original FASTA file %s" %f[6:]
	with open (f, "r") as fasta_file:
		## parse input to get a clean file_id
		sequences = fasta_file.read()
		fasta_file.close()
		print "Done."

		print "Starting BLASTN Search against %s" %f[6:]
		print "This may take awhile, and requires a presistant Internet connection..."
		## Blast using the FASTA set, and assing the XML string result to result_handle
		result_handle = NCBIWWW.qblast("blastn", "nt", sequences, expect=8.0)
		print "Done."
		print "Writing results to a BLAST.XML file..."
		## write the XML string to a file for later parsing
		blast_file_id = "".join(["blast-temp/", f[6:-6], ".xml"])
		save_result = open(blast_file_id, "w")
		save_result.write(result_handle.read())
		save_result.close()
		result_handle.close()
		print "Done."


	print "BLAST search successfully completed.\n"

print "Process Complete."

