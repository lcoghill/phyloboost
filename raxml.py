from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO
from Bio import AlignIO
import glob
import os

##get list of all new fasta files

print "\n\nGetting a list of alignment files..."
alignment_files = glob.glob("alignments/*.phylip") # get a list of all fasta files in /alignment
file_count = len(alignment_files)
print "%s files successfully found.\n" %file_count
print "-"*50+"\n"


## need to convert the files to phylip format so that RAxML can read them cleanly
working_dir = os.path.dirname(os.path.realpath(__file__)) + "/trees/"

for f in alignment_files:
	print "Building Trees for %s" % f
	raxml_cline = RaxmlCommandline(sequences=f, model="GTRCAT", name=f[12:-7], working_dir=working_dir)
	stdout, stderr = raxml_cline()
	print "Done.\n"
