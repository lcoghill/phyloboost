from Bio import SeqIO
import glob

## get list of all new fasta files

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("output/expanded-fasta/*.fasta")  # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successful found.\n" % file_count

bad_char_list = []  # just used for bug tracking.

for f in fasta_files:  # iterate through all the files in the directory
    print "\nChecking FASTA file %s for non-nucleotide characters..." % f
    clean_seq_list = []  # build a new list of good sequences
    bad_seq_count = 0  # keep track of how many sequences are discarded
    tracker = 0  # track if a sequence is bad or not.
    for seq_record in SeqIO.parse(f, "fasta"):  # for each sequence

        dna = seq_record.seq  # assign it to a string

        for char in dna:  # for each character in the sequence
            if char not in (
                    'A', 'T', 'C', 'G', 'N', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'V', 'H', 'D', 'X', 'a', 't', 'c',
                    'g',
                    'n', 'u', 'r', 'y', 'k', 'm', 's', 'w', 'b', 'v', 'h', 'd',
                    'x'):  # if its not an accepted base or ambiguity code
                bad_char_list.append(char)  # add it to the bad list for bug tracking
                tracker = 1  # set the bad / good tracker to bad

        if tracker == 1:  # if bad sequence
            bad_seq_count += 1  # increment bad counter
        elif tracker == 0:  # if good sequence
            clean_seq_list.append(seq_record)  # add the record to the good list
        tracker = 0  # after done iterating through a sequence, reset tracker for the next sequence
    print "Complete."
    if bad_seq_count == 0:
        print "There were no non-nucleotide characters in FASTA file %s." % f
    else:
        print
        'There were %s non-nucleotide containing entries found in the FASTA file. Those entries have been removed.' \
        % bad_seq_count

    SeqIO.write(clean_seq_list, f,
                "fasta")  # after all iterations, write the good seq record list to the same file we started with

print "\nFASTA clean complete."
