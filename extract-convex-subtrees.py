import ivy
from ivy import treegraph as tg

def taxid_count (tree):

    records = tree.split(",")
    tis = [r.split('_ti')[1].replace(")", "") for r in records]

    return tis

def color_vertices(taxonomy, treegraph, tid):
    """
    tid: NCBI taxon id
    
    Color the vertices of `treegraph` that are members of taxon `tid`
    """
    nxt, bck = taxonomy.hindex[taxonomy.taxid_vertex[tid]]

    colored = treegraph.new_vertex_property('bool')
    # `colored` is a boolean vertex property map that will flag those
    # vertices in the unrooted treegraph that are in convex subtrees
    # corresponding to taxon `tid`
    
    seen = set()
    lvs = set()
    for v in treegraph.vertices():
        if v.out_degree() == 1: # leaf
            seen.add(v)
            taxv = taxonomy.taxid_vertex[treegraph.vertex_taxid[v]]
            if taxonomy.incertae_sedis[taxv]:
                p = taxv.in_neighbours().next()
                pn, pb = taxonomy.hindex[p]
                if nxt >= pn and bck <= pb:
                    colored[v] = 1
                    lvs.add(v)
            else:
                n, b = taxonomy.hindex[taxv]
                if n >= nxt and b <= bck:
                    colored[v] = 1
                    lvs.add(v)

    def gather():
        s = set()
        for v in lvs:
            for n in v.out_neighbours():
                if not colored[n]: s.add(n)
        return s

    def check(x):
        i = 0
        for y in x.out_neighbours():
            if not colored[y]: i += 1
        return i

    verts = gather()
    while 1:
        for x in verts:
            if check(x) == 1:
                lvs.add(x)
                colored[x] = 1
        v = gather()
        if v == verts: break
        verts = v

    c = tg.defaultdict(list)
    # `c` is a mapping of positive integer values to lists of colored
    # vertices in treegraph. The integers are counts of adjacent
    # vertices that are not colored. So if `taxv` corresponds to a
    # convex subgraph of `treegraph`, `c` should be {1: [x]}, where x
    # is the vertex point of attachment
    for v in treegraph.vertices():
        if colored[v]:
            i = check(v)
            if i: c[i].append(v)

    return colored, c

def proc(g, line, merged, probfile, outfile):
    pbtree, s = line.split()
    print 'processing', pbtree
    r = ivy.newick.parse(s) # the root node of the tree of interest
    lvs = r.leaves()
    rps = [] # rootpaths of leaf nodes, where each rootpath is a list
             # of taxids from leaf to root
    leaf_tid_counts = tg.Counter()
    try:
        for lf in lvs:
            # assign/compute attributes of leaves
            w = lf.label.split('_')
            lf.gi = int(w[-2][2:])
            lf.taxid = int(w[-1][2:])
            leaf_tid_counts[lf.taxid] += 1
            if lf.taxid not in g.taxid_vertex and lf.taxid in merged:
                lf.taxid = merged[lf.taxid]
            ## lf.taxv = g.taxid_vertex[lf.taxid]
            taxv = g.taxid_vertex[lf.taxid]
            lf.taxid_next, lf.taxid_back = g.hindex[taxv]
            lf.taxid_rootpath = tg.taxid_rootpath(g, lf.taxid)
            for i, x in enumerate(lf.taxid_rootpath):
                if x not in g.taxid_vertex and x in merged:
                    lf.taxid_rootpath[i] = merged[x]
            rps.append(lf.taxid_rootpath)
    except:
        print '!!! problem assigning leaf taxids'
        probfile.write('%s\n' % pbtree)
        #return []

    r.mrca = tg.rootpath_mrca(rps) # taxid of mrca of all tree's leaves

    taxids = set()
    for rp in rps:
        # trim rootpaths: make them terminate with mrca
        while 1:
            if rp[-1] == r.mrca: break
            else: rp.pop()
        assert rp
        taxids.update(rp)

    # create a taxonomy (sub)graph of only those taxids represented in r
    ## taxidsubg = tg.taxid_subgraph(g, taxids)
    taxidsubg = tg.taxid_new_subgraph(g, taxids)
    taxidsubg.vfilt = taxidsubg.new_vertex_property('bool')

    ## r.taxv = taxidsubg.taxid_vertex[r.mrca]

    # no need to check for convexity for singleton tip taxa
    for x in [ taxidsubg.taxid_vertex[lf.taxid] for lf in lvs
               if leaf_tid_counts[lf.taxid]==1 ]:
        taxidsubg.vfilt[x] = 0
    
    # an undirected graph having the same topology as r, used for
    # checking convexity of taxa
    treegraph = tg.gt.Graph(directed=False)
    treegraph.mrca = r.mrca
    print 'mrca:', g.taxid_name(r.mrca)
    treegraph.vertex_taxid = tg.get_or_create_vp(treegraph, 'taxid', 'int')
    ## treegraph.vertex_taxv = tg.get_or_create_vp(treegraph, 'taxv', 'int')
    v2lf = {}
    N = len(r)
    verts = treegraph.add_vertex(N)
    for n in r: # for each node in r
        # store its treegraph vertex
        n.v = verts.next()
        if not n.children:
            treegraph.vertex_taxid[n.v] = n.taxid
            ## treegraph.vertex_taxv[n.v] = int(n.taxv)
            v2lf[n.v] = n
        if n.parent:
            treegraph.add_edge(n.parent.v, n.v)

    treegraph_leaves = [ x for x in treegraph.vertices() if x.out_degree()==1 ]
    convex = {} # for storing the convex subgraphs
    def traverse(taxv):
        """
        `taxv` is a vertex in the taxonomy graph. This function checks whether
        it is convex in `treegraph`; if yes, stores the info in
        `convex`; if no, it recursively checks descendants of `taxv` for
        convexity
        """
        tid = taxidsubg.vertex_taxid[taxv]
        print 'checking', tid, taxidsubg.vertex_name[taxv]
        p, c = color_vertices(g, treegraph, tid)
        if len(c)==1 and len(c[1])==1: # taxv/tid is convex
            print '...success'
            rv = c[1][0] # rv is the root of the convex subtree
            treegraph.set_vertex_filter(p)
            ## lvs = [ x for x in treegraph.vertices() if x.out_degree()==1 ]
            lvs = [ x for x in treegraph_leaves if p[x] ]
            if len(lvs) > 2:
                # we are only interested in convex subgraphs having
                # more than 2 leaves
                rootpaths = []
                for lf in lvs:
                    ti = treegraph.vertex_taxid[lf]
                    tv = taxidsubg.taxid_vertex[ti]
                    if not taxidsubg.incertae_sedis[tv]:
                        rootpaths.append(tg.taxid_rootpath(taxidsubg, ti))
                if rootpaths:
                    mrca = tg.rootpath_mrca(rootpaths)
                    print 'traverse: mrca', mrca
                    ancv = [taxidsubg.taxid_vertex[mrca]]
                    while ancv[-1] != taxv:
                        # STRANGE EDGE CASES HERE
                        try: ancv.append(ancv[-1].in_neighbours().next())
                        except StopIteration: pass

                    k = '.'.join([ str(taxidsubg.vertex_taxid[x])
                                   for x in ancv ])
                    convex[k] = (rv, p)
            treegraph.set_vertex_filter(None)
        else:
            treegraph.set_vertex_filter(None)
            for n in taxv.out_neighbours():
                traverse(n)

    for v in taxidsubg.root.out_neighbours(): traverse(v)
    ## print 'done'

    def make_newick(root, seen):
        children = [ x for x in root.out_neighbours() if x not in seen ]
        if children:
            seen.update(children)
            s = '(%s)' % ','.join(
                [ make_newick(c, seen) for c in children ]
                )
        else:
            s = v2lf[root].label.replace(',','').replace('(','').replace(')','')
        return s
        
    newicks = []
    for k, (root, p) in convex.items():
        treegraph.set_vertex_filter(p)
        s = make_newick(root, set([root]))
        treegraph.set_vertex_filter(None)
        names = ','.join([ g.taxid_name(int(x)) for x in k.split('.') ])
        tis = taxid_count(s)
        print set(tis)
        print len(set(tis))
        if len(set(tis)) > 2:
            outfile.write('%s\t%s\t%s\t%s;\n' % (pbtree, k, names, s))
            print 'wrote subtree:', names

    for n in r.postiter():
        n.parent = None; del n.children

if __name__ == "__main__":
    merged = {}
    with open('ncbi/merged.dmp') as f:
        for line in f:
            v = line.split()
            merged[int(v[0])] = int(v[2])

    g = tg.load_taxonomy_graph('ncbi/ncbi.xml.gz')

    probfile = open('readable.problem_subtrees','w')
    outfile = open('readable.convex_subtrees','w')
    with open('trees.out') as f:
        for line in f:
            proc(g, line, merged, probfile, outfile)
    outfile.close()
    probfile.close()
