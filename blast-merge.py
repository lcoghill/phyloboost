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
gi_stats = [];
seq_stats = [];

def fasta_process(file_name, blast_file):
	global new_seqs
	global orig_seqs
	global gi_stats
	global seq_stats
	sequences = []
	orig_sequence_gis = []
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
	print "Complete."
	new_list = []
	gi_final = []

	for g in gi_list:
		
		string = "".join(g)
		string = string.encode("ascii")
		new_list = string.split("|")
		gi_final.append(new_list.pop(0))


	## Use the list of gi values to retreive the appropriate genbank FASTA data for those records
	print "Getting original sequences from FASTA file %s" %file_id[2:]
	with open (f, "r") as fasta_file:
		## parse input to get a clean file_id		
		orig_sequences = fasta_file.readlines()
		orig_seqs = orig_seqs + len(orig_sequences)
		fasta_file.close()
		orig_sequence_gis = [s for s in orig_sequences if ">gi" in s]
		print "Complete."

	orig_gi_count = 0
	for o in orig_sequence_gis:
		orig_sequence_gis[count] = o[4:].strip("\n")
		gi_final.append(orig_sequence_gis[count])
		orig_gi_count = orig_gi_count + 1
		
	print "Grabbing genbank data for all GIs..."
	new_fasta_file_name = "".join(["fasta/expanded-fasta/",file_id,".FASTA"])
	new_fasta_file = open(new_fasta_file_name, "w+")
	
	gi_final = list(set(gi_final)) ## remove any duplicate GI values in the list.

	for gi in gi_final:
		genbank_raw = library.genbank.fetchseq(gi)
		sequences.append("".join([">gi|", gi, "\n"]))
		sequences.append("".join([str(genbank_raw.seq), "\n"]))	

	
	largest_seqs = []
	smallest_seqs = []
	largest_seq_length = len(max(sequences)) - 1
	smallest_seq_length = largest_seq_length
	total_seq_count = len(sequences)/2
	for line in sequences:
		if len(line) == len(max(sequences)):
			current_largest_id = sequences[int(sequences.index(line))-1]
			largest_seqs.append(current_largest_id[4:-2])
		if line[0] != ">" and (len(line)-1) < largest_seq_length and (len(line)-1) < smallest_seq_length:
			smallest_seq_length = len(line) - 1
			current_small_id = sequences[int(sequences.index(line))-1] #If a match is found, append to the smallest_seqs list after reducing the current index by one. Getting the GI instead of the sequence.
			smallest_seqs.append(current_small_id[4:-2])
		new_fasta_file.write("%s" % line)

	seq_stats.append(largest_seq_length)
	seq_stats.append(largest_seqs)
	seq_stats.append(smallest_seq_length)
	seq_stats.append(smallest_seqs)
	seq_stats.append(total_seq_count)	
	print "Complete."
	new_seqs = len(gi_final);
	gi_stats = gi_final;
	return (new_seqs, orig_seqs, gi_stats, seq_stats);







### Program Code ###

fasta_files = glob.glob("fasta/*.FASTA") # get a list of all fasta files in /fasta
blast_files = glob.glob("blast-results/*.xml") # get a list of all fasta files in /fasta

count = 0;

stats_file = open("output_stats.txt", "a+")

for f in fasta_files:
	
	## process each file by calling fasta_process
	fasta_process(f, blast_files[count])
	print "FASTA file %s successfully merged." % f
	count = count + 1
	total_new_seqs = new_seqs + total_new_seqs;
	total_orig_seqs = orig_seqs + total_orig_seqs;
	ti_ci = f.split("_")
	ti_ci[0] = ti_ci[0][8:]
	ti_ci[1] = ti_ci[1][2:-6]
	gi_string = ",".join(gi_stats)
	stats_string = "".join(["Filename: ",f,"\nTi: ",ti_ci[0],"\nCi: ",ti_ci[1],"\nStarting Sequences: ",str((orig_seqs/2)),"\nNew Sequences: ",str(new_seqs),"\nTotal Count: ",str(seq_stats[4]),"\nLongest Sequence(s): ",(",".join(seq_stats[1])),"\nLength: ",str(seq_stats[0]),"\nShortest Sequence(s): ",(",".join(seq_stats[3])),"\nLength: ",str(seq_stats[2]),"\nList of new GIs: ","\n",gi_string,"\n","#" * 30,"\n"])
	stats_file.write(stats_string);
	orig_seqs = 0
	seq_stats = []
print "Complete."
print "Adding additional records to FASTA complete."
print "New FASTA files can be found in /fasta/expanded-fasta directory.\n\n\n"
print "############### SNAPSHOT STATS ###############\n"
print "Total Files Processed: %s" %len(fasta_files)
print "Total Original Sequences: %s" % (total_orig_seqs/2)
print "Total Additional Sequences Included: %s" % total_new_seqs
print "\n############### SNAPSHOT STATS ###############"
stats_file.close();
