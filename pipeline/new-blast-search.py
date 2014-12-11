from multiprocessing import Pool
from Bio import SeqIO
import numpy as np
import subprocess
import glob
import os



def median_seq_length(f, sequences) :
	
	lengths = []
	sequence_ids = []
	handle = open(f, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		lengths.append(len(record.seq))
		sequence_ids.append(record.id[3:])
	handle.close()

	median_seq_len = np.median(lengths)
	min_len = int(median_seq_len - (median_seq_len * length_filter))
	max_len = int(median_seq_len + (median_seq_len * length_filter))
	
	return min_len, max_len, sequence_ids

def blast_search(orig_f, e_value, blast_db, length_filter, blast_files, threads) :
	
	min_len, max_len, sequence_ids = median_seq_length (f, length_filter)
	out_file = blast_files + f.split("/ti")[1]
	subprocess.call(["blastn", "-query", str(orig_f), "-db", str(blast_db), "-evalue", 
		str(e_value), "-outfmt", "6", "-out", str(out_file),
		"-num_threads", str(threads), "-use_index", "true"])
	## function here to purge duplicate GIS and hits < 100% identity and make sure there are only the correct number
	clean_results(sequence_ids, orig_num_seqs)
	## call blastdbcmd here to get sequences for the appropriate number of hits
	subprocess.call(["blastdbcmd", "", "", "", ""])
	## function here to combine to FASTA files into single expanded-fasta file
	build_fasta()
	## remove blast results hits file
	os.remove(blast_hits_file)

def clean_results(sequence_ids, blast_hits_file) :
	
	hits = []
	handle = open(blast_hits_file, 'r')
	for line in handle :
		hits.append(line.strip())
	handle.close()

	good_hits = []
	for hit in hits :
		if hit[1] == 100 and hit[0] not in sequence_ids :
			good_hits.append(hit)

	handle = open(blast_hits_file, 'w')
	if len(good_hits) + len(sequence_ids) <= 5000 :
		for hit in good_hits :
			handle.write(hit+"\n")
			handle.close()
	else :
		num_hits = 5000 - len(sequence_ids)
		short_hits = good_hits[:num_hits]
		for hit in short_hits :
			handle.write(hit+"\n")
			handle.close()
	
def build_fasta(orig_f, blast_fasta) :
	## open blast fasta
	## open original fasta
	## combine two
	## write over blast fasta



def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return ([e for e in t if e != None] for t in itertools.izip_longest(*args))




blast_db = "/home/lcoghill/Dev/phyloboost/pipeline/phyloboost_1.5_blast/phyloboost_1.5"   # local copy of blast nt database
cluster_dir = "/home/lcoghill/Dev/phyloboost/pipeline/clusters/"   # location where fasta files generated from get-clusters.py are stored. Must have .fasta suffix on each file.
blast_files = "/home/lcoghill/Dev/phyloboost/pipeline/expanded-fasta/"    # location where blast result files will be saved.
e_value = 10                                  # e-value threshold for blast search
length_filter = .50                           # length in bp to filter blast results at. (hits longer than this value will be discarded)
blast_threads = 8


## Get a count of the number of fasta files in the appropriate directory
print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("".join([cluster_dir,"*.fas"]))
file_count = float(len(fasta_files))
print "%s files successfully found.\n" %file_count

for f in fasta_files :
	print "BLASTING file %s..." %f
	blast_search(f, e_value, blast_db, length_filter, blast_files, blast_threads)
### read in a fasta file
### calc the median seq length
### perform blast search
### keep records that meet criteria
### write new fasta file with new records
