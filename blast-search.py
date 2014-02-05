import os
import re
import sets
import glob
import logging
import MySQLdb
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successfully found.\n" %file_count

## WARN IF THERE ARE A LOT OF FASTA FILES
if file_count > 2:
	print "WARNING:\n"
	print "There are %s FASTA files to be processed." %file_count
	print "This process can take quite a long time and requires a presistant Internet connection...\n\n"

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
		result_handle = NCBIWWW.qblast("blastn", "nt", sequences, expect=8.0)
		print "Done."
		print "Writing results to a BLAST.XML file..."
		## write the XML string to a file for later parsing
		blast_file_id = "".join(["blast-temp/", f[6:-6], ".xml"])
		save_result = open(blast_file_id, "w")
		save_result.write(result_handle.read())
		save_result.close()
		result_handle.close()
		print "Step %s / %s complete." % (count + 1, file_count)
		count = count + 1


	print "BLAST search successfully completed.\n"

print "Process Complete."

