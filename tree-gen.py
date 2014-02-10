from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO
from Bio import AlignIO
import glob

##get list of all new fasta files

print "\n\nGetting a list of alignment files..."
alignment_files = glob.glob("alignments/*.align") # get a list of all fasta files in /alignment
file_count = len(alignment_files)
print "%s files successful found.\n" %file_count


## need to convert the files to phylip format so that RAxML can read them cleanly

for f in alignment_files:

	input_handle = open(f, "rU")
	print f[11:-6]
	alignments = AlignIO.read(input_handle, "clustal")
	AlignIO.write(alignments, f[11:-6], "phylip-relaxed")
	input_handle.close()


## due to the -n limitation of no "/" characters, must and pass it as an argument?

#in_file = "alignments/testfile.phylip"

#raxml_cline = RaxmlCommandline(sequences=in_file, model="GTRCAT", name="test-tree")
#stdout, stderr = raxml_cline()

#input_handle.close()

