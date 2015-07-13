import ivy
from ivy import treegraph as tg
from ivy.tree import Node
from collections import Counter

def pull_line(filename, i):
    "utility function for testing purposes - pull a single line from a file"
    with open(filename) as f:
        for i in range(i):
            line = f.next()
    return line

def lkup(x):
    "utility function for translating taxids into names"
    if isinstance(x, int):
        return g.taxid_name(x)
    elif isinstance(x, str):
        return g.taxid_name(int(x.split('_')[1][2:]))
    elif isinstance(x, Node):
        if x.label:
            return lkup(x.label)
        return x
    # not int, str (node label), or node; maybe a collection
    try:
        return [ lkup(e) for e in x ]
    except TypeError:
        pass
    return x

def name2taxid(name):
    return [ g.vertex_taxid[v] for v in g.vertices()
             if g.vertex_name[v] == name ]

def trace_conflicting_nodes(r, taxid):
    """
    given r (the root of a mapped tree) and a conflicting taxid, return
    the node(s) that collectively contradict the monophyly of 'taxid'
    """
    lvs = r.conflicts[taxid]
    mrca = r.con2mrca[taxid]
    c = Counter()
    for lf in lvs:
        for anc in lf.rootpath(end=mrca):
            if anc is not mrca:
                c[anc] += 1
    return [ anc for anc, i in c.items() if i < len(anc.children) ]

def categorize_taxa_and_nodes(g, r):
    """

    g: taxonomic hierarchy (from ivy.treegraph.load_taxonomy_graph)

    r: root node (ivy.tree.Node) of a tree

    The information content of a phylogenetic tree is necessarily
    expressed in terms of names (taxa) arranged in a hierarchy
    (taxonomy). A rooted tree relating sampled members of taxa
    (individual exemplar organisms) 'represents' a subset of higher taxa
    in the hierarchy.

    Definition: a taxon is 'represented' in a tree if it appears in the
    sub-hierarchy that connects all the taxa assigned to leaf nodes to
    their MRCA in the taxonomy.

    Here we map (align) a tree to a given taxonomy and categorize the
    taxa represented in it. Each taxon falls in 1 (and only 1) of the
    following:
    
     1. monophyletic - 2 or more immediate descendants of the taxon
        (or members of a terminal taxon) are represented in the tree
        and form an exclusive clade. In practice, these are the taxa
        in position 0 of the node.taxids array of internal nodes.
    
     2. singleton - only 1 immediate descendant of the taxon (or
        member of a terminal taxon) is represented. In practice,
        these are taxa that appear in position 1 or higher in the
        node.taxids array of internal nodes, or position 0 of leaf
        nodes.
    
        3. contradicted - taxa with 2 or more immediate descendants (or
        members of a terminal taxon) represented in the tree, and those
        descendants/members do not form an exclusive clade; but more
        stringently, it is not possible to make the taxon monophyletic
        simply by resolving (soft) polytomies, because 1 or more
        strongly supported nodes break up the descendants/members into
        separate clades. A contradicted taxon maps to >= 1 such node in
        the tree.
    
     4. 'softly' contradicted - taxa with 2 or more immediate
        descendants (or members of a terminal taxon) represented in
        the tree that are not monophyletic, but could be made
        monophyletic by resolving a single polytomy.
    
    Only categories 1 and 3 are interesting (informative)
    
    
    Nodes: each node in the tree falls in only 1 category:
    
     1. leaf node
    
     2. internal node that represents a monophyletic taxon and zero
        or more ancestral singleton taxa
    
     3. internal node that represents a union set of taxa
        (node.stem_cdef) that does not contradict any taxon's
        monophyly (i.e., it resolves relationships among descendants
        of a taxon)
    
     4. internal node that represents a union set of taxa
        (node.stem_cdef) that contradicts the taxonomic hierarchy in
        some way -- a contradicted taxon's leaf members will trace
        through >= 1 of these nodes on their way to their mrca in
        the tree
    
     5. root node
    
    Only categories 2, 3, and 4 are interesting (informative)

    This function maps the tree to the taxonomy and returns dictionaries
    mapping taxa (taxids) to nodes. Summary statistics such as numbers
    of taxa or nodes in each category can be computed by counting the
    keys or values in these dictionaries.
    
    """
    # map the tree to taxonomy
    r.ladderize()
    ivy.tree.index(r)
    for n in r.leaves():
        gistr, tistr, accstr = n.label.split('_')
        n.snode_id = int(gistr[2:]) # gi
        n.taxid = int(tistr[2:]) # ti
    tg.map_stree(g, r)

    mono_taxid2node = {} # map monophyletic taxid -> node
    singleton_taxid2node = {} # map singleton taxid -> node
    for x in r:
        if x.taxids:
            if x.children:
                mono_taxid2node[x.taxids[0]] = x
            else:
                for taxid in x.taxids[1:]:
                    singleton_taxid2node[taxid] = x

    contradicted_taxid2nodes = {} # map contradicted taxid -> [node(s)]
    softly_contradicted_taxid2nodes = {} # map softly contradicted taxid -> node
    for x in r.conflicts.keys():
        v = trace_conflicting_nodes(r, x)
        if v:
            contradicted_taxid2nodes[x] = v
        else:
            softly_contradicted_taxid2nodes[x] = r.con2mrca[x]
    
    return (mono_taxid2node, singleton_taxid2node, contradicted_taxid2nodes,
            softly_contradicted_taxid2nodes)
    
def compute_summary_stats(tree_outfile, taxa_outfile, treeid, clusterid, r,
                          mono_taxid2node, singleton_taxid2node,
                          contradicted_taxid2nodes,
                          softly_contradicted_taxid2nodes):
    """

    Write summary stats to outfile.

    Stats are designed to help answer questions like:

      Which trees are the most informative? the least informative?
      where 'informative' is calculated in terms of taxa supported and
      contradicted as monophyletic.

      Which taxa are most frequently supported? contradicted?

      Which trees are the most/least informative in terms of 'new'
      (unnamed) clades that are compatible/incompatible with taxonomy?

    """

    ntax_mono = len(mono_taxid2node)
    ntax_singleton = len(singleton_taxid2node)
    ntax_contradicted = len(contradicted_taxid2nodes)
    ntax_softly = len(softly_contradicted_taxid2nodes)

    leaves = r.leaves()
    nleaves = len(leaves)
    
    # number of internal nodes, excluding root -- i.e. potentially
    # informative nodes)
    nnodes = len([ x for x in r if x.children and x.parent ])

    # 'conflict' nodes (traced by members of contradicted taxa to
    # their mrca in the tree): raw number, and weighted by number of
    # taxa contradicted
    conflict_nodes = Counter()
    for v in contradicted_taxid2nodes.values():
        for x in v:
            conflict_nodes[x] += 1
    nnodes_conflicting = len(conflict_nodes)
    nnodes_conflicting_weighted = sum(conflict_nodes.values())

    # nodes that resolve relationships without contradicting taxonomy
    nnodes_resolving = nnodes - ntax_mono - nnodes_conflicting

    # number of taxa for which the tree is potentially informative
    ntax_total = ntax_mono + ntax_contradicted + ntax_softly
    
    if ntax_total > 0 and nnodes > 0 :
	    # taxon monophyly stats:
	    for taxid in mono_taxid2node:
	        name = g.taxid_name(taxid)
	        row = [taxid, name, treeid, 1]
	        line = '\t'.join(map(str, row))
	        taxa_outfile.write('{}\n'.format(line))

	    for taxid in contradicted_taxid2nodes:
	        name = g.taxid_name(taxid)
	        row = [taxid, name, treeid, 0]
	        line = '\t'.join(map(str, row))
	        taxa_outfile.write('{}\n'.format(line))

	    nnodes_float = float(nnodes)
	    ntax_total_float = float(ntax_total)

	    # tree stats
	    row = [
	        treeid,
	        clusterid,
	        r.taxid, # taxid of tree's mrca
	        g.taxid_name(r.taxid),
	        ntax_total, # number of taxa about which the tree is potentially
	                    # informative
	        ntax_mono, # no. taxa supported == no. nodes supporting monophyly
	        ntax_contradicted,
	        ntax_softly,
	        ntax_singleton,
	        nnodes, # number of potentially informative nodes
	        nnodes_conflicting, # no. nodes contradicting taxa
	        nnodes_conflicting_weighted, # weighted number of above
	        nnodes_resolving,
	        # proportion of nodes supporting monophyly
	        ntax_mono/nnodes_float,
	        # proportion of nodes contradicting monophyly
	        nnodes_conflicting/nnodes_float,
	        # proportion of 'resolving' nodes (neither supporting nor
	        # contradicting taxonomy)
	        nnodes_resolving/nnodes_float,       
	        # proportion of 'positive' nodes
	        (ntax_mono + nnodes_resolving)/nnodes_float,
	        # proportion of taxa supported as monophyletic
	        ntax_mono/ntax_total_float,
	        # proportion of taxa contradicted as monophyletic
	        ntax_contradicted/ntax_total_float,
	        # how skeletonized is the tree?
	        (ntax_total + ntax_singleton)/nnodes_float
	        ]
	    line = '\t'.join(map(str, row))
	    tree_outfile.write('{}\n'.format(line))

if __name__ == '__main__':
    #import vis
    
    print "Reading taxonomy graph..."
    g = tg.load_taxonomy_graph('ncbi.gt.gz')
    
    trees_file = 'readable.convex_subtrees.out'

    ## i = 99
    ## line = pull_line(trees_file, i)
    ## v = line.split('\t')
    ## clusterid = v[0]
    ## stree = v[-1]
    ## r = ivy.tree.read(stree)
    ## vis.view(g, r)

    tree_outfile = open('tree_stats.csv','w')
    tree_headers = (
        'treeid clusterid taxid name '
        'ntax ntax_mono ntax_contradiced ntax_softly ntax_singleton '
        'nnodes nnodes_conflicting nnodes_conflicting_weighted nnodes_resolving '
        'nnodes_prop_mono nnodes_prop_conflicting nnodes_prop_resolving '
        'nnodes_prop_positive ntax_prop_mono ntax_prop_contradicted skeletal'
        ).replace(' ', '\t')
    tree_outfile.write('{}\n'.format(tree_headers))
    
    taxa_outfile = open('taxa_stats.csv','w')
    taxa_headers = 'taxid name treeid mono'.replace(' ','\t')
    taxa_outfile.write('{}\n'.format(taxa_headers))

    with open(trees_file) as f:
        for treeid in range(10):
            line = f.next()
            v = line.split('\t')
            clusterid = v[0]
            stree = v[-1]
            r = ivy.tree.read(stree)
            (mono_taxid2node, singleton_taxid2node, contradicted_taxid2nodes,
             softly_contradicted_taxid2nodes) = categorize_taxa_and_nodes(g,r)
            
           compute_summary_stats(
                tree_outfile, taxa_outfile, treeid, clusterid, r,
                mono_taxid2node, singleton_taxid2node,
                contradicted_taxid2nodes, softly_contradicted_taxid2nodes)
