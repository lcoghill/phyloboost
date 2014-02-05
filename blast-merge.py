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
### Functions ###

new_seqs = 0;
total_new_seqs = 0;
total_orig_seqs = 0;
orig_seqs = 0;

def fasta_process(file_name, blast_file):
	global new_seqs
	global orig_seqs

	with open (file_name, "r") as fasta_file:
		## parse input to get a clean file_id
		file_id = file_name.split("/")
		file_id = file_id[1].split(".")
		file_id = file_id[0]
		

	### Parses the XML Blast return document, and returns a clean list of only gi (genebank ids) for the sequences of interest
	print "Parsing BLAST Results file for GI values..."
	for record in NCBIXML.parse(open(blast_file)) :
	    #We want to ignore any queries with no search results:
	    gi_list = []
	    if record.alignments :
	        for align in record.alignments :
	            for hsp in align.hsps :
	                if hsp.expect < 0.04 and align.length <= 7500: #if the result has an expected value less than 0.04 (strong match) and is less than 7500bp in length (max length from Phylota) keep it.
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
	print "Getting original sequences from FASTA file %s" %file_id[2:]
	with open (f, "r") as fasta_file:
		## parse input to get a clean file_id		
		sequences = fasta_file.readlines()
		orig_seqs = orig_seqs + len(sequences)
		fasta_file.close()
		print "Done."
	print "Grabbing genbank data for new GIs..."
	new_fasta_file_name = "".join(["fasta/expanded-fasta/",file_id,".FASTA"])
	new_fasta_file = open(new_fasta_file_name, "w+")
	
	for line in sequences:
 		 new_fasta_file.write("%s" % line)		
		
	for gi in gi_final:
		genbank_raw = library.genbank.fetchseq(gi)
		fasta_entry = str("".join([">gi|", gi, "\n", str(genbank_raw.seq), "\n"]))		
		new_fasta_file.write(fasta_entry)

	new_fasta_file.close()
	print "Done."
	new_seqs = new_seqs + len(gi_final);
	return (new_seqs, orig_seqs);

### Program Code ###

fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
blast_files = glob.glob("blast-temp/*.xml") # get a list of all fasta files in /fasta

count = 0;


for f in fasta_files:
	
	## process each file by calling fasta_process
	fasta_process(f, blast_files[count])
	print "FASTA file %s successfully merged." % f
	count = count + 1
	total_new_seqs = new_seqs + total_new_seqs;
	total_orig_seqs = orig_seqs;

print "Complete."
print "Adding additional records to FASTA complete."
print "New FASTA files can be found in /fasta/expanded-fasta directory.\n\n\n"
print "############### SNAPSHOT STATS ###############\n"
print "Total Files Processed: %s" %len(fasta_files)
print "Total Original Sequences: %s" % (total_orig_seqs/2)
print "Total Additional Sequences Included: %s" % total_new_seqs
print "\n############### SNAPSHOT STATS ###############"
