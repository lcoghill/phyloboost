from ivy import treegraph as tg
from collections import Counter
import graph_tool.all as gt
from Bio import SeqIO
import numpy as np
import glob



def mad(data) :
    ## calulate the median absolute deviation for each distance from the mcp
    ## median of the set comprising the absolute values of the differences between 
    ## the median and each data point.
    return np.median(np.abs(data - np.median(data)))

def flag_outliers(data, keys, cutoff) :
    ## flag any sequence as an outlier that is more than 'cutoff' standard deviations
    ## away from median distance of the rest of the points.
    
    outliers = []
    if len(set(data)) > 1 :
        mad_res = mad(data)
        ## adjusts for the MAD=0 problem.
        ## if MAD=0 then all sequences are kept, and the file will be flagged for manual review.
        if mad_res > 0 :
            mad_data = np.abs(data - np.median(data)) / mad_res
        else :
            mad_data = []    
        
        for i, m in enumerate(mad_data) :
            if m > cutoff :
                outliers.append(keys[i])

    return outliers

def fetch_gis(cluster) :
    gis_with_tis = {}
    sequences = list(SeqIO.parse(open(cluster), 'fasta'))
    for s in sequences :
        rec_id = s.id.split("_")
        gi = int(rec_id[0][2:])
        ti = int(rec_id[1][2:])
        gis_with_tis[gi] = ti

    return gis_with_tis

def write_seqs(cluster, outliers, out_cluster_dir) :
    sequences = list(SeqIO.parse(open(cluster), 'fasta'))
    file_name = out_cluster_dir + cluster.split("/")[-1]
    good_seqs = []
    outlier_seqs = []
    out_handle = open(file_name, 'a')
    for s in sequences :
        gi = int(s.id.split("_ti")[0][2:])
        
        if gi not in outliers :
            good_seqs.append(s)
        else :
            outlier_seqs.append(s)

    count = SeqIO.write(good_seqs, out_handle, 'fasta')
    out_handle.close()

    return count, outlier_seqs

def log_outliers(outlier_seqs, log_file) :
    out_handle = open(log_file, 'a')
    for s in outlier_seqs :
        rec_id = s.id.split("_")
        gi = rec_id[0][2:]
        ti = rec_id[1][2:]
        out_handle.write(",".join([gi,ti,str(mcp[0][0]),str(s.seq), "\n"]))

    out_handle.close()



## needed globals variables
in_cluster_dir = ''
out_cluster_dir = ''
ncbi_graph_file = 'ncbi.gt.gz'
log_file = 'outlier_log.csv'
## this is the number of standard deviations away from the distribtuion median 
## a sequence must be before it will be flagged and removed.
cutoff = 3.0  

## set for testing 
#gis_with_tis = {1: 7994, 2: 178766, 3: 643502, 4: 4362, 5: 270330, 6: 930229}

clusters = glob.glob(in_cluster_dir + "*.fas")
print "Found %i clusters." %len(clusters)
print "Loading NCBI taxonomy graph..."
g = tg.load_taxonomy_graph(ncbi_graph_file)

status = 1
total = len(clusters)
good_seq_count = 0
bad_seq_count = 0
out_handle = open(log_file, 'a')
out_handle.write(",".join(["GI", "TI", "MCP", "Sequence", "\n"]))
out_handle.close()

if __name__ == '__main__' :

    for c in clusters :
        print "Processing cluster %s, cluster %i / %i..." %(c, status, total)
        gis_with_tis = fetch_gis(c)
        cas = []
        check = {}
        ## get the root path of each gi in the dict
        ## fetch the parent of each gi
        ## reset the graph to directed to use tg.taxid_rootpath
        g.set_directed(True)
        for key, val in gis_with_tis.items() :
        ## find the parent of each TI in the graph
        ## find the most common parent for the entire cluster
            if g.taxid_vertex[val] :
                cas.append(tg.taxid_rootpath(g, val)[1])
                

        ## count occurences to find the most common parent
        count = Counter(cas)
        mcp = count.most_common(1)
        

        ## have to remove the direction of the graph in order to get a distance
        ## from one tip to another. Otherwise in directed graphs, shortest_distance
        ## can only be calculated in the same direction as the graph
        g.set_directed(False)

        ## calculate the shortest distance from each tip, to the most-common-parent
        distances = {}
        for key, val in gis_with_tis.items() :
            source = g.taxid_vertex[val]
            target = g.taxid_vertex[mcp[0][0]]
            dist = gt.shortest_distance(g, source, target)
            distances[key] = dist

        ## book-keeping 
        dst = []
        keys = []
        for key, val in distances.items() :
            keys.append(key)
            dst.append(val)

        ## flag outliers that are more than "cutoff" SDs away from group median
        outliers = flag_outliers(dst, keys, cutoff)
        
        ## write non-outlier sequences to new fasta file, and log outlier sequences for further review
        good_count, outlier_seqs = write_seqs(c, outliers, out_cluster_dir)
        print "A total of %i sequences passed; a total of %i suspected outlier sequence(s) were removed." %(good_count, len(outlier_seqs))

        ## some logging variables and functions
        good_seq_count = good_seq_count + good_count
        bad_seq_count = bad_seq_count + len(outliers)
        ## write the outlier sequences
        log_outliers(outlier_seqs, log_file)

        status += 1

print "-"*100
print "Check complete."
print "\n%i clusters checked." %len(clusters)
print "%i sequences seem to be properly identified." %good_seq_count
print "%i possible mis-identified sequence(s) removed from clusters.\n" %bad_seq_count
print "-"*100