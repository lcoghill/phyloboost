import ivy



in_file = 'trees.out'
out_file = '70_support_collapsed_trees.out'
collapse_value = 70.00



in_handle = open(in_file, 'r')
out_handle = open(out_file, 'a')
tree_count = sum(1 for line in open(in_file))
prog_count = 1
bad_trees = []
for line in in_handle :
    record = line.strip("\n").split("\t")
    print "Collapsing nodes with support under %s percent for file %s, tree %s / %s....." %(collapse_value, record[0], prog_count, tree_count)
    tree = ivy.tree.read(record[-1])
    leaves = tree.leaves()
    for x in tree :
        if not x.isroot and x not in leaves and x.label is None :
            x.label = 0.00
    orig_node_cnt = len( [ x for x in tree if x not in leaves and not x.isroot ])
    low_support_nodes = [ x for x in tree if x not in leaves and not x.isroot and float(x.label) <= collapse_value ]
    [ x.collapse() for x in reversed(low_support_nodes) ]    
    newick = tree.write()
    
    try :
        new_tree = ivy.tree.read(newick)
        new_leaves = new_tree.leaves()
        check_set = [ x.label for x in new_tree if x not in new_leaves and not x.isroot and x.label is None ]
        if len(check_set) == 0 :    
            check = True
    except :
        bad_trees.append(line.strip())
        check = False
    
    if check :
        out_str = "\t".join([record[0], record[1], record[2], newick])
        out_handle.write(out_str + "\n")
    
    prog_count += 1

in_handle.close()
out_handle.close()

if not bad_trees :
    print "All trees seem to be able to be parsed by IVY."
else :
    print "There are trees with incorrect newick formatting."
    print [ x for x in bad_trees ]
