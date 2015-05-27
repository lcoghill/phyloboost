from ivy import treegraph as tg
from Bio import SeqIO
import numpy as np
import glob



def mad(data) :
    ## calulate the median absolute deviation for each distance from the mcp
    ## median of the set comprising the absolute values of the differences between 
    ## the median and each data point.
    return np.median(np.abs(data - np.median(data)))

def flag_outliers(distances, cutoff) :
    ## flag any sequence as an outlier that is more than 'cutoff' standard deviations
    ## away from median distance of the rest of the points.
    keys = []
    vals = []
    for key, val in distances.items() :
        keys.append(key)
        vals.append(val)

    outliers = []
    ## attempts to adjusts for the MAD=0 problem.
    if len(set(vals)) > 1 :
        mad_res = mad(vals)
        if mad_res > 0 :
            mad_data = np.abs(vals - np.median(vals)) / mad_res
        else :
            mad_data = []    
        
        for i, m in enumerate(mad_data) :
            if m > cutoff :
                outliers.append(keys[i])

    return outliers

def fetch_gis(cluster, bad_tis) :
    gis_with_tis = {}
    sequences = list(SeqIO.parse(open(cluster), 'fasta'))
    for s in sequences :
        rec_id = s.id.split("_")
        gi = int(rec_id[0][2:])
        ti = int(rec_id[1][2:])
        
        ## checking that all tis are in taxonomy graph, if not, flagging and logging for manual
        ## review and removal from clusters
        if not g.taxid_name(ti) :
            bad_tis.append(ti)
        else :    
            gis_with_tis[gi] = ti

    return gis_with_tis, bad_tis

def calc_distances(gis_with_tis, g) :
    
    ## gathers the rootpath for each ti in a cluster
    ## then looks for the intersection for each rootpath with all other rootpaths in a cluster
    ## calculates the number of steps to each TI in the result of the intersection
    ## keeps the shortest distance
    
    gis = []
    tis = []
    for key, val in gis_with_tis.items() :
        tis.append(val)
        gis.append(key)

    distances = {}
    rootpaths = []
    for t in tis :
        rootpaths.append(tg.taxid_rootpath(g, t))

    for i, r in enumerate(rootpaths) :
        ccas = []
        matches = []
        for r2 in rootpaths :
            if r != r2 :
                paths = []
                paths.append(r)
                paths.append(r2)
                tcca = tg.rootpath_mrca(paths)
                matches.append(tcca)
       
        if len(set(matches)) > 1 :
            dist_set = []
            for m in set(matches) :
                dist_set.append(r.index(m))
            dist = min(dist_set)
            cca = r[dist_set.index(dist)]
        else :
           cca = matches[0]
           dist = r.index(cca)
        
        distances[gis[i]] = dist

    return distances

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

def log_outliers(c, outlier_seqs, log_file) :
    out_handle = open(log_file, 'a')
    clust = c.split("/")[-1]
    for s in outlier_seqs :
        rec_id = s.id.split("_")
        gi = rec_id[0][2:]
        ti = rec_id[1][2:]
        out_handle.write(",".join([clust, gi, ti, str(s.seq), "\n"]))

    out_handle.close()



## needed globals variables
in_cluster_dir = ''
out_cluster_dir = ''
ncbi_graph_file = 'ncbi.gt.gz'
log_file = 'outlier_log.csv'
## this is the number of standard deviations away from the distribtuion median 
## a sequence must be before it will be flagged and removed.
cutoff = 5.0  

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
out_handle.write(",".join(["Cluster", "GI", "TI", "Sequence", "\n"]))
out_handle.close()



if __name__ == '__main__' :

    for c in clusters :
        bad_tis = []
        print "Processing cluster %s, cluster %i / %i..." %(c, status, total)
        gis_with_tis, bad_tis = fetch_gis(c, bad_tis)
        cas = []
        check = {}
        ## get the root path of each gi in the dict
        ## fetch the parent of each gi
        ## reset the graph to directed to use tg.taxid_rootpath
        g.set_directed(True)

        ## calculate the distance for each TI in a cluster to the closest shared ancestor with any
        ## other TI in a cluster.
        distances = calc_distances(gis_with_tis, g)
        
        ## flag outliers that are more than "cutoff" SDs away from group median
        outliers = flag_outliers(distances, cutoff)
        
        ## write non-outlier sequences to new fasta file, and log outlier sequences for further review
        good_count, outlier_seqs = write_seqs(c, outliers, out_cluster_dir)
        print "A total of %i sequences passed; a total of %i suspected outlier sequence(s) were removed." %(good_count, len(outlier_seqs))

        ## some logging variables and functions
        good_seq_count = good_seq_count + good_count
        bad_seq_count = bad_seq_count + len(outliers)
        ## write the outlier sequences
        log_outliers(c, outlier_seqs, log_file)

        status += 1
        ## dump missing ti list to file to check / remove from clusters
        bad_tis_h = open('bad_tis.txt', 'a')
        for t in bad_tis :
            bad_tis_h.write(str(t) + "\n")
        bad_tis_h.close()


print "-"*100
print "Check complete."
print "\n%i clusters checked." %len(clusters)
print "%i sequences seem to be properly identified." %good_seq_count
print "%i possible mis-identified sequence(s) removed from clusters.\n" %bad_seq_count
print "-"*100