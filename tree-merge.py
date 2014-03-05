<BS>from StringIO import StringIO
from Bio import SeqIO
from Bio import AlignIO
import glob

#get a list of all "best" tree files in the directory

##get list of all new fasta files

print "\n\nGetting a list of tree files..."
fasta_files = glob.glob("trees/expanded-fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successful found.\n" %file_count


## read in each tree




#append the ci and ti in phylota format
## write each set to a new file
