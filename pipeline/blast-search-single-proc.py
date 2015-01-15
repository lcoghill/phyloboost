from Bio import SeqIO
import numpy as np
import subprocess
import itertools
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

def blast_search(f, blast_files, e_value, length_filter, blast_db, cores) :

	file_name = f.split("/")[-1].replace(".fas", "")
	min_len, max_len, sequence_ids = median_seq_length (f, length_filter)
	out_file = blast_files + f.split("/")[-1].replace(".fas", ".tsv")

	subprocess.call(["megablast", "-i", str(f), "-d", str(blast_db), "-e", str(e_value),
		"-a", str(cores), "-D", "3", "-I", "T", "-f", "T", "-p", "98.0", "-o", str(out_file)])
		
	clean_results(sequence_ids, out_file)
	subprocess.call(["blastdbcmd", "-db", str(blast_db), "-dbtype", "nucl", "-entry_batch", out_file, 
		"-out", out_file.replace(".tsv", ".fas"), "-outfmt", "%f"])
		
	build_fasta(f, out_file, max_len, min_len)
	os.remove(out_file)
		
def clean_results(seq_accs, blast_hits_file) :
	
	hits = []
	hit_accs = []
	handle = open(blast_hits_file, 'r')
	for line in handle :
		if not line.startswith("#") :
			hit_rec = line.split()
			hit_acc = hit_rec[1]
			if hit_acc not in hit_accs and hit_acc not in seq_accs :
				hit_accs.append(hit_acc)
				hits.append(line.strip())
	
	handle.close()

	handle = open(blast_hits_file, 'w')
	if len(hits) + len(seq_accs) <= 5000 :
		for hit in hits :
			rec = hit.split()
			handle.write(rec[1]+"\n")
		handle.close()
	else :
		num_hits = 5000 - len(seq_accs)
		short_hits = hits[:num_hits]
		for hit in short_hits :
			rec = hit.split()
			handle.write(rec[1]+"\n")
		handle.close()
	
def build_fasta(orig_f, out_file, max_len, min_len) :
	blast_log = open("blast.log", "a")
	expanded_fasta = []
	blast_handle = open(out_file.replace(".tsv", ".fas"), "r")
	orig_handle = open(orig_f, "r")
	blast_fasta = list(SeqIO.parse(blast_handle, "fasta"))
	added_recs = len(blast_fasta)
	len_filtered_recs = 0
	for n in blast_fasta :
		if len(n.seq) <= max_len and len(n.seq) >= min_len :
			n.id = "gi|" + n.id.split("|")[1]
			n.name = ""
			n.description = ""
			expanded_fasta.append(n)
		else :
			len_filtered_recs += 1

        orig_fasta = list(SeqIO.parse(orig_f, "fasta"))
	blast_log.write(orig_f + "\t|\t" + str(len(orig_fasta)) + "\t|\t" + str(added_recs - len_filtered_recs) + "\n")
	blast_log.close()
	for o in orig_fasta :
		expanded_fasta.append(o)

	output_handle = open(out_file.replace(".tsv", ".fas"), "w")
	SeqIO.write(expanded_fasta, output_handle, "fasta")
	output_handle.close()

def div_list(fasta_files, cores) :
	cores = int(cores)
	return [ fasta_files[i::cores] for i in xrange(cores) ]

def find_new_files(fasta_files, cluster_dir) :

	print "Searching blast results for already expanded fasta files..."
        handle = open('blast.log', 'r')
        done_names = []
        new_files = []
        for line in handle :
            name = line.strip().split()[0].replace(cluster_dir, "")
            done_names.append(name)
	for f in fasta_files :
		name = f.replace(cluster_dir, "")
		if name not in done_names :
			new_files.append(name)

        handle.close()
	return new_files




blast_db = "phyloboost_1.5_blast/phyloboost_1.5"   # local copy of blast nt database
cluster_dir = ""   # location where fasta files generated from get-clusters.py are stored. Must have .fasta suffix on each file.
blast_output_dir = ""    # location where blast result files will be saved.
e_value = 10                                  # e-value threshold for blast search
length_filter = .50                           # length in bp to filter blast results at. (hits longer than this value will be discarded)
cores = 4

## Get a count of the number of fasta files in the appropriate directory
print "\n\nGetting a list of FASTA files..."
all_fasta_files = glob.glob("".join([cluster_dir,"*.fas"]))
file_count = float(len(all_fasta_files))
print "%s files successfully found.\n" %file_count

if os.path.exists(os.getcwd() + '/blast.log') :
    print "Found previous blast log file."
    fasta_files = find_new_files(all_fasta_files, cluster_dir)
else :
    fasta_files = all_fasta_files


count = 1
if len(fasta_files) > 0 :
    print "There appears to be %s fasta files without blast results logged." %len(fasta_files)
    for f in fasta_files :
        print "Blasting file %s / %s..." %(count, file_count)
    	blast_search(f, blast_output_dir, e_value, length_filter, blast_db, cores)
        count += 1
print "All Done."
