from Bio import SeqIO
import glob



clusters_dir = '' # location of cluster files from vsearch
pic_dir = '' # directory to save the phylogenetically informative clusters


clusters = glob.glob(clusters_dir)

for c in clusters :
    sequences = SeqIO.parse(open(c), 'fasta')
    tis = []
    for s in sequences :
        ti = s.id.split("_")[1][2:]
        tis.append(ti)
    if len(set(tis)) > 3 :
        out_file = pic_dir + c.split("/")[-1]
        out_handle = open(out_file, 'a')
        SeqIO.write(sequences, out_handle, 'fasta')
        out_handle.close()
