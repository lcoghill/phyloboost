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

## Need to add feature to iterate through all files in the alignment directory.
## Also need to investigate the ability to control where the RAxML output is saved.
input_handle = open("alignments/testfile.align", "rU")
output_handle = open("alignments/testfile.align", "w+) 
alignments = AlignIO.read(input_handle, "clustal")
AlignIO.write(alignments, output_handle, "phylip-relaxed")


## due to the -n limitation of no "/" characters, must and pass it as an argument?

in_file = "alignments/testfile.phylip"

raxml_cline = RaxmlCommandline(sequences=in_file, model="GTRCAT", name="test-tree")
stdout, stderr = raxml_cline()

input_handle.close()

