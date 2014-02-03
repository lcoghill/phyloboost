from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO

input_handle = open("alignments/testfile.align.FASTA", "rU")
output_handle = open("testfile.nex", "w")

for record in SeqIO.parse(input_handle, "fasta") :
	sequences = SeqIO.parse(input_handle, "fasta")
	print record.id
	SeqIO.write(sequences, output_handle, "nexus")

output_handle.close()
handle.close()

#raxml_cline = RaxmlCommandline(sequences="testfile.nex", model="PROTCATWAG", name="interlaced2")
#stdout, stderr = raxml_cline()
