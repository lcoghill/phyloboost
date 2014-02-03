from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO
from Bio import AlignIO

align = AlignIO.read("alignments/testfile.align", "clustal")
print align

#raxml_cline = RaxmlCommandline(sequences="testfile.nex", model="PROTCATWAG", name="interlaced2")
#stdout, stderr = raxml_cline()
