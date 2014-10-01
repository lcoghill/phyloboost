from StringIO import StringIO
from Bio import AlignIO
from Bio import SeqIO
import numpy as np
import glob


def get_median(f):
    lengths = []
    rec_count = 0
    for seq_record in SeqIO.parse(f, "fasta"):
        lengths.append(len(seq_record.seq))
        rec_count = rec_count + 1
    median = np.median(lengths)
    return (median, rec_count)

def correct_taxid(ti, ti_dict) :

    if ti in ti_dict :
        new_ti = ti_dict[ti]
        return new_ti
    else:
        return ti




fasta_dir = "" # full path to the directory of fasta files
ti_error_file = 'merged.dmp' # usually merged.dmp from NCBI
## builds a dict of ncbi deprecated ti values and their replacements
ti_handle = open(ti_error_file, 'r')
ti_dict = {}
for line in ti_handle :
    ti_dict[line.split("|")[0].replace("\t", "")] = line.split("|")[1].replace("\t", '')


print "*"*100
print "This script modifies a directory of fasta files in the following ways:"
print " - Removes sequences that contain non-IUPAC codes." 
print " - Removes sequences that are longer or shorter than a given threshold of file median."
print " - If you wish to skip length filtering, set the percentage value to 0."
print "*"*100

## Get list of all new fasta files
print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob(fasta_dir + "*.fasta") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successful found.\n" %file_count

percent = float(raw_input("Please enter a percentage threshold as a decimal (ie: 0.50):"))

bad_char_list = [] # just used for bug tracking.

for f in fasta_files: # iterate through all the files in the directory
    print "\nChecking FASTA file %s for non-nucleotide characters..." %f
    clean_seq_list = [] # build a new list of good sequences
    bad_seq_count = 0 # keep track of how many sequences are discarded
    tracker = 0 # track if a sequence is bad or not.
    for seq_record in SeqIO.parse(f, "fasta"): # for each sequence
        ti = seq_record.id.split("_")[1][2:]
        gi = seq_record.id.split("_")[0][2:]
        ti = correct_taxid(ti, ti_dict)
        seq_record.id = "".join(["gi",gi,"_ti",ti])
        seq_record.description = ""
        dna = seq_record.seq.lower() # assign it to a string
        for char in dna: # for each character in the sequence
            if char not in ('a','t','c','g','n','u','r','y','k','m','s','w','b','v','h','d','x'): # if it is not an accepted base or ambiguity code
                bad_char_list.append(char) # add it to the bad list for bug tracking
                tracker = 1 # set the bad / good tracker to bad
    
        if tracker == 1: # if bad sequence
            bad_seq_count = bad_seq_count + 1 #increment bad counter
        elif tracker == 0: # if good sequence
            clean_seq_list.append(seq_record) #add the record to the good list
        tracker = 0 # after done iterating through a sequence, reset tracker for the next sequence
    print "Complete."
    if bad_seq_count == 0:
        print "There were no non-nucleotide characters in FASTA file %s." % f
    else:
        print "There were %s non-nucleotide containing entries found in the FASTA file. Those entries have been removed." % bad_seq_count
    
    SeqIO.write(clean_seq_list, f, "fasta") # after done with all iterations, write the good seq record list to the same file we started with

if percent > 0 :
    for f in fasta_files: # iterate through all the files in the directory
        print "\nChecking FASTA file %s for sequences of unusual lengths..." %f
        count = 0
        new_records = []
        median, rec_count = get_median(f)
        min_len = median-(percent*median)
        max_len = median+(percent*median)
        print "Removing sequences that are shorter or longer than %s%% (%s, %s) of the median length (%s bp)..." %(percent*100, min_len, max_len, median)
        for seq_record in SeqIO.parse(f, "fasta"):
            if len(seq_record.seq) > max_len or len(seq_record.seq) < min_len:
                count = count + 1.
            else:
                new_records.append(seq_record)

        SeqIO.write(new_records, f, "fasta")
        print "%s sequences removed out of %s total sequences from file %s due to unsual lengths." %(count, rec_count, f)
        

print "\nFASTA clean complete."
