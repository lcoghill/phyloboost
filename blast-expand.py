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


### Functions

def fasta_process(file_name):
	with open (file_name, "r") as fasta_file:
		## parse input to get a clean file_id
		file_id = file_name.split("/")
		file_id = file_id[1].split(".")
		file_id = file_id[0]

	### Parses the XML Blast return document, and returns a clean list of only gi (genebank ids) for the sequences of interest
	print "Parsing BLAST Results file for GI values..."
	for record in NCBIXML.parse(open(blast_file_id)) :
	    #We want to ignore any queries with no search results:
	    gi_list = []
	    if record.alignments :
	        for align in record.alignments :
	            for hsp in align.hsps :
	                if hsp.expect < 0.04: #if the result has an e-value less than .04, keep the result and put it in a new fasta sequence file with all appropriate information ****NEED TO FIGURE OUT HOW TO LIMIT LENGTH ****
	                	gi_list.append(align.hit_id.split("gi|"))
	print "Done."
	new_list = []
	gi_final = []
	for g in gi_list:
		
		string = "".join(g)
		string = string.encode("ascii")
		new_list = string.split("|")
		gi_final.append(new_list.pop(0))


	### Use the list of gi values to retreive the appropriate genbank FASTA data for those records
	print "Grabbing genbank data for new GIs..."
	new_fasta_file_name = "".join(["fasta/expanded-fasta/",file_id,".FASTA"])
	new_fasta_file = open(new_fasta_file_name, "w+")
	new_fasta_file.write(sequences)
		
	for gi in gi_final:
		genbank_raw = library.genbank.fetchseq(gi)
		fasta_entry = str(">gi|" + gi + "|" + (genbank_raw.description) + "\n" + (genbank_raw.seq) + "\n")  ##replace with a join for speed and to make more concise		
		new_fasta_file.write(fasta_entry)

	new_fasta_file.close()
	print "Done."

### Program Code

fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta

for f in fasta_files:
	## process each file by calling fasta_process
	fasta_process(f)
	print "FASTA file %s successfully expanded." % f

print "Cleaning up undeed files..."
## remove uneeded BLAST XML Result Files

print "Done."

print "Adding additional records to FASTA complete."
print "New FASTA files can be found in /fasta/expanded-fasta directory."

