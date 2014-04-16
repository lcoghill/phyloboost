from Bio import SeqIO
from Bio.Phylo.Applications import RaxmlCommandline
from StringIO import StringIO
from Bio import AlignIO
import glob
import os
import re

## locate trees already generated
def done_trees(trees.directory):
    print "Getting a list of trees already built..."
    tree_files = glob.glob(trees.directory) # get a list of all files in trees.directory
    complete_trees = []
    for t in tree_files:
        name = t[28:]
        complete_trees.append(name) 
    return complete_trees



## write tree files to single output tree file
def write_trees(tree.outfile, trees.directory):

    tree_files = glob.glob("".join([trees.directory,"RAxML_bestTree.*"]))
    tree_file = open(tree.outfile, 'w+')
    for f in tree_files:
        print "-"*50
        print "Writing trees to file $s" %tree.outfile
        newick_string = open(f)
        record = "".join([name,"\t",newick_string,"\n"])
        tree_file.write(record) 
        print "Tree %s, successfully written." %f



## delete all tree files when done?
def delete_trees(keep_all_trees, trees.directory):

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
tree.outfile = "trees.out"                          # final, combined output file name
trees.directory = "/trees"                          # directory where tree files should be saved
keep_all_trees = 1                                  # flag whether to keep all tree output files when the process is done. (0 = Delete Everything, 1 = Keep Everything)
alignment.files = "/alignments"                     # location where all alignment files are stored. (Must be in Phylip format with .phy extension)
num.threads = 1                                     # number of threads to use in RAxML. Don't use more than physical cores in on the computer.
num.replicates = 1                                  # number of replicates for RAxML
evo.model = "GTRCAT"                                # model of evolution used for RAxML

###





complete_trees = done_trees(trees.directory) # get a list of the trees already built, by comparing files in working dir vs those in alignment.files
alignments = glob.glob("".join([alignment.files,".phy"])) # get a list of all phylip files in /alignment
alignment_names = [] #names of alignment files without directory or suffix
print "\n\nGetting a list of alignment files..."
for a in alignments:
    match = re.match(r"^.*\/(.*).phy.*$",a)
    name = match.group(1)
    alignment_names.append(name)

s = set(complete_trees)
files_to_use = [x for x in alignment_names if x not in s] ## keep the files that are in the alignment directory, but not in the complete_trees list
print "%s files new files found.\n" % len(files_to_use)
print "-"*50+"\n"


## build the trees using RAxML
count = 0 ## status counter
for f in files_to_use:
    infile = "".join([alignments,f,'.phy'])
    print "Builing tree for %s..." %f
    raxml_cline = RaxmlCommandline(sequences=infile, model=evo.model, name=f, working_dir=trees.directory, threads=num.threads, num_replicates=num.replicates) ## raxml wrapper call. Assumes the use of the pthreads version of raxml 8.0 or above
    stdout, stderr = raxml_cline()
    print "Tree %s / %s Complete." %(count, len(files_to_use))
    count += 1 

write_trees(tree.outfile) # write all "best" trees to a single file in newick format
delete_trees(keep_all_trees, trees.directory) # delete or keep all tree files
tree_file.close()
print "-"*50
print "Tree file generation complete."
