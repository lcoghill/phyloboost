import MySQLdb
import os
import logging
import library.genbank
logging.basicConfig(level=logging.INFO)
from Bio import Entrez, SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
import re
import sets



### A few import static variables
library.genbank.email = 'me@my.address.com'    


with open ("fasta/ti2763_ci2.FASTA", "r") as fasta_file:
	sequences = fasta_file.read()
	fasta_file.close()
	#print sequences

#     result_handle = NCBIWWW.qblast("blastn", "nt", sequences)
#     save_result = open("blast_result.xml", "w")
#     save_result.write(result_handle.read())
#     save_result.close()
#     result_handle.close()

### Parses the XML Blast return document, and returns a clean list of only gi (genebank ids) for the sequences of interest

for record in NCBIXML.parse(open("blast_result.xml")) :
    #We want to ignore any queries with no search results:
    gi_list = []
    if record.alignments :
        for align in record.alignments :
            for hsp in align.hsps :
                if hsp.expect < 0.04: #if the result has an e-value less than .04, keep the result and put it in a new fasta sequence file with all appropriate information ****NEED TO FIGURE OUT HOW TO LIMIT LENGTH ****
                	gi_list.append(align.hit_id.split("gi|"))

new_list = []
gi_final = []
for g in gi_list:
	
	string = "".join(g)
	string = string.encode("ascii")
	new_list = string.split("|")
	gi_final.append(new_list.pop(0))


### Use the list of gi values to retreive the appropriate genbank FASTA data for those records
new_fasta_file = open("fasta/expanded-fasta/ti2763_ci2.FASTA", "w+")
new_fasta_file.write(sequences)
	
for gi in gi_final:
	genbank_raw = library.genbank.fetchseq(gi)
	fasta_entry = str(">gi|" + gi + "|" + (genbank_raw.description) + "\n" + (genbank_raw.seq) + "\n")  ##replace with a join for speed and to make more concise		
	new_fasta_file.write(fasta_entry)

new_fasta_file.close()