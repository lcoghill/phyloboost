from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO
from Bio import AlignIO
import glob
import os
import re

## locate trees already generated
def done_trees(trees_directory):
    print "Getting a list of trees already built..."
    tree_files = glob.glob(trees_directory) # get a list of all files in trees.directory
    complete_trees = []
    for t in tree_files:
        name = t[28:]
        complete_trees.append(name) 
    return complete_trees



## write tree files to single output tree file
def write_trees(tree_outfile, trees_directory):

    tree_files = glob.glob("".join([trees_directory,"RAxML_bestTree.*"]))
    tree_file = open(tree_outfile, 'w+')
    for f in tree_files:
        print "-"*50
        print "Writing trees to file $s" %tree_outfile
        newick_string = open(f)
        record = "".join([name,"\t",newick_string,"\n"])
        tree_file.write(record) 
        print "Tree %s, successfully written." %f



## delete all tree files when done?
def delete_trees(keep_all_trees, trees_directory):

    if keep_all_trees == 0:
        # delete all files in trees directory here
        treefile_list = glob.glob("".join([trees.directory,"*.*"]))
        for f in treefile_list:
            os.remove(f)
    elif keep_all_trees == 1:
        print "Keeping all tree files."
    else:
        print "Inproper keep_all_trees flag. Keeping all tree files just to be safe."





### Some important variables                        ###
tree_outfile = "trees.out"                          # final, combined output file name
trees_directory = ""                          # directory where tree files should be saved
keep_all_trees = 1                                  # flag whether to keep all tree output files when the process is done. (0 = Delete Everything, 1 = Keep Everything)
alignment_dir = ""                     # location where all alignment files are stored. (Must be in Phylip format with .phy extension)
num_threads = 4                                     # number of threads to use in RAxML. Don't use more than physical cores in on the computer.
num_replicates = 100                                  # number of replicates for RAxML
evo_model = "GTRCAT"                                # model of evolution used for RAxML
p_seed = 12345 # parsimony seed for RAxML
bs_seed = 12345 # bootstrap seed for RAxML
algorithm = 'a' # bootstrap algorithm for RAxML

###



complete_trees = done_trees(trees_directory) # get a list of the trees already built, by comparing files in working dir vs those in alignment.files
alignments = glob.glob(alignment_dir + '/*.phy') # get a list of all phylip files in /alignment_dir
alignment_names = [] #names of alignment files without directory or suffix
print "\n\nGetting a list of alignment files..."
for a in alignments:
    name = re.search('alignments/(.+?).phy', a).group(1)
    alignment_names.append(name)

s = set(complete_trees)
files_to_use = [x for x in alignment_names if x not in s] ## keep the files that are in the alignment directory, but not in the complete_trees list
print "%s files new files found.\n" % len(files_to_use)
print "-"*50+"\n"


## build the trees using RAxML
count = 0 ## status counter
for f in files_to_use:
    infile = "".join([alignment_dir, "/", f,'.phy'])
    print infile
    print "Builing tree for %s..." %f
    raxml_cline = RaxmlCommandline(sequences=infile, model=evo_model, name=f, working_dir=trees_directory, threads=num_threads, num_replicates=num_replicates, algorithm=algorithm, rapid_bootstrap_seed=bs_seed, parsimony_seed=p_seed) ## raxml wrapper call. Assumes the use of the pthreads version of raxml 8.0 or above
    stdout, stderr = raxml_cline()
    print "Tree %s / %s Complete." %(count, len(files_to_use))
    count += 1 

write_trees(tree_outfile, trees_directory) # write all "best" trees to a single file in newick format
delete_trees(keep_all_trees, trees_directory) # delete or keep all tree files
print "-"*50
print "Tree file generation complete."
