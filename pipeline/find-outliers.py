from ivy import treegraph as tg
from collections import Counter
from glob import glob
import numpy as np

ncbi_graph_file = 'ncbi.gt.gz' #taxonomy graph file
cluster_dir = '' #directory of vsearch produced clusters
cutoff = 3 # number of mean absolute devations from median

g = tg.load_taxonomy_graph(ncbi_graph_file)
for clustfile in glob(cluster_dir):
    headers = [ x[1:-1].split('_') for x in open(clustfile) if x[0]=='>' ]
    gi2ti = dict([ (int(a[2:]), int(b[2:])) for a,b in headers ])
    tis = sorted(set(gi2ti.values()))
    rootpaths = [ tg.taxid_rootpath(g, ti) for ti in tis ]
    ## print 'mrca:', g.taxid_name(tg.rootpath_mrca(rootpaths))
    counts = Counter()
    for rp in rootpaths:
        for ti in rp:
            counts[ti] += 1

    def f(rp):
        'steps to most recent common ancestor of any other ti in the cluster'
        for i, ti in enumerate(rp):
            if counts[ti]>1:
                break
        return i

    steps = [ f(rp) for rp in rootpaths ]
    median = np.median(steps)
    absdev = [ abs(x-median) for x in steps ]
    mad = np.mean(absdev)

    res = [ ti for x, ti in zip(absdev, tis) if x > cutoff ]
    if res:
        print clustfile, res
