from Bio import Phylo
import csv
import numpy
import dendropy
from dendropy import treecalc

def to_distance_matrix(tree):
    """Create a distance matrix (NumPy array) from clades/branches in tree.
 
    A cell (i,j) in the array is the length of the branch between allclades[i]
    and allclades[j], if a branch exists, otherwise infinity.
 
    Returns a tuple of (allclades, distance_matrix) where allclades is a list of
    clades and distance_matrix is a NumPy 2D array.
    """
    allclades = list(tree.find_clades(order='level'))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    distmat = numpy.repeat(numpy.inf, len(allclades)**2)
    distmat.shape = (len(allclades), len(allclades))
    for parent in tree.find_clades(terminal=False, order='level'):
        for child in parent.clades:
            if child.branch_length:
                distmat[lookup[parent], lookup[child]] = child.branch_length
    if not tree.rooted:
        distmat += distmat.transpose
    return (allclades, numpy.matrix(distmat))




tree_id = []
trees = []
csv.field_size_limit(300000)

import csv
with open('test-tree.txt', 'rb') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        tree_id.append(row[0])
        trees.append(row[1])

trees = Phylo.parse(trees, 'newick')
