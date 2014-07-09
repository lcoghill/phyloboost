from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import csv
import os


def save_blast(results_file):
    save = raw_input('Would you like to save the blast results file? (Y/N):')
    if save == 'N' or save == 'n':
        os.remove(results_file)
        print "Results not saved. Exiting now."
    else:
        print "File saved as %s." %results_file

def add_to_fasta(query_file, results_file, closest_cluster, original_fasta):
    new_file = closest_cluster + ".fasta.new"
    raw_fasta_file = original_fasta + new_file.strip(".new")
    
    query_handle = open(query_file, 'rU')
    query_record = list(SeqIO.parse(query_handle, "fasta"))
    original_handle = open(raw_fasta_file, 'rU')
    records = list(SeqIO.parse(original_handle, "fasta"))
    output_handle = open(new_file, 'w')
    new_records = []

    for r in records:
        new_records.append(r) 
    for r in query_record:
        new_records.append(r)

    SeqIO.write(new_records, output_handle, "fasta")
    output_handle.close()
    output_handle.close()
    query_handle.close()
    original_handle.close() 
    
    print "New fasta saved as %s." %new_file
    
    save_blast(results_file)


blast_db = 'phylota_blastdb_v1.0'
query_file = 'test.fasta'
e = 0.001
results_file = 'results.csv'
original_fasta_dir = 'fasta/'

print "\n"

print "Blasting target file %s against database %s..." %(query_file, blast_db)

blastn_cline = NcbiblastnCommandline(query=query_file, db=blast_db, evalue=e, outfmt=6, max_target_seqs=1, out=results_file)
stdout, stderr = blastn_cline()

with open(results_file, 'rb') as csvfile:
    csv_file = csv.reader(csvfile, delimiter='\t')
    for row in csv_file:
        print row
        closest_cluster = row[1].split("_cluster_")[1]
        source_id = row[0]
        percent_match = row[2]
        score = row[11]

print "Process Complete.\n"
print "#"*75
print "Source ID: %s" %source_id
print "Closest Match: %s" % closest_cluster
print "Percent Match: %s" % percent_match
print "E-value: %s" %e
print "Score: %s" %score
print "#"*75
print "\n"

choice = raw_input("Would you like to add the sequence " + source_id + " to the fasta file for cluster "+ closest_cluster+"? (Y/N):")

if choice == 'Y' or choice == 'y':
    add_to_fasta(query_file, results_file, closest_cluster, original_fasta_dir)
elif choice == 'N' or choice == 'n':
    save_blast(results_file)      
else:
    print "Invalid choice, blast results are saved to %s." %results_file
