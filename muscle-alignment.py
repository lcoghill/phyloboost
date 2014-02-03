from Bio.Align.Applications import MuscleCommandline
from StringIO import StringIO
from Bio import AlignIO
import glob

##get list of all new fasta files

print "\n\nGetting a list of FASTA files..."
fasta_files = glob.glob("fasta/expanded-fasta/*.FASTA") # get a list of all fasta files in /fasta
file_count = len(fasta_files)
print "%s files successful found.\n" %file_count


### iterate through each file doing the following:

for f in fasta_files: 
	print "Aligning FASTA file %s" % f
	file_name = f[21:-6]
	alignment = MuscleCommandline(input=f, out="".join(["alignments/",file_name,".align"]), clwstrict=True)
	stdout, stderr = alignment()
	print "Success."

