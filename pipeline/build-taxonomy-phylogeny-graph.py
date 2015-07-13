import ivy, requests
from ivy import treegraph as tg
from collections import defaultdict
import graph_tool.all as gt



taxonomy_graph = '' # name / location of taxonomy graph file 
out_file_name = '' # name of output graph file


g = tg.load_taxonomy_graph(taxonomy_graph)
strees = []
i = 1
with open('convex_subtrees.out') as f:
    for line in f:
        tree = line.split("\t")[-1]
        print tree
        r = ivy.tree.read(tree)
        ivy.tree.index(r)
        for n in r:
            if n.isleaf:
                v = n.label.split('_')
                n.snode_id = int(v[0][2:])
                n.taxid = int(v[1][2:]) if (len(v)>1 and
                                        v[1] and v[1] != 'None') else None
            else:
                n.snode_id = int(n.id)
        r.stree = i
        strees.append(r)
        i += 1

stree2color = {}
for i, r in enumerate(strees):
    stree2color[r.stree] = tg.color20[i % 20]

taxids = set()
for r in strees:
    tg.map_stree(g, r)
    for lf in r.leaves(): taxids.update(lf.taxid_rootpath)

root_taxa = set([ r.taxid for r in strees ])
if len(root_taxa) > 1:
    rps = [ tg.taxid_rootpath(g, x) for x in root_taxa ]
    mrca_taxid = tg.rootpath_mrca(rps)
    for x in rps:
        taxids.update(x[:x.index(mrca_taxid)+1])

taxg = tg.taxid_new_subgraph(g, taxids)

verts = taxg.new_vertex_property('bool')
edges = taxg.new_edge_property('bool')
for r in strees:
    tg.merge_stree(taxg, r, r.stree, verts, edges)
 
root = taxg.root

gv = tg.graph_view(taxg, vfilt=verts, efilt=edges)

for x in gv.vertices():
    if x.in_degree()==0 and int(x)!=int(taxg.root):
        v = taxg.vertex(int(x))
        while 1:
            try:
                e = v.in_edges().next()
            except StopIteration:
                break
            edges[e] = 1
            p = e.source()
            if verts[p]: break
            else:
                verts[p] = 1
                v = p

gv.save(out_file_name)
