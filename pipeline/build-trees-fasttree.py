from Bio import AlignIO
import subprocess
import dendropy
import glob
import os
import re



def cal_gap_prop(alignment_dir, file_name) :
    aln = AlignIO.read(alignment_dir + file_name + ".fas", 'fasta')
    gap_count = 0
    for a in aln :
        gap_count = gap_count + a.seq.count("-")
    total_count = len(aln) * len(aln[1].seq)
    ad = round(1.0 - float(gap_count) / float(total_count), 2)
    print ad

    return ad

## locate trees already generated
def done_trees(trees_directory):
    print "Getting a list of trees already built..."
    tree_files = glob.glob(trees_directory) # get a list of all files in trees.directory
    complete_trees = []
    for t in tree_files:
        name = t[28:]
        complete_trees.append(name) 
    return complete_trees

## fix support values in newick
def fix_newick(newick) :
    itree = ivy.tree.read(newick)
    leaves = itree.leaves()
    for x in itree :
        if not x.isroot and x not in leaves and x.label is not None :
            new_label = float(x.label)
            x.label = int(new_label*100)
    fixed_newick = itree.write()
    
    return fixed_newick

## write tree files to single output tree file
def write_trees(file_name, tree_outfile, trees_directory, tree_type):


    tree_file = "".join([trees_directory, file_name])
    tl = dendropy.Tree.get_from_path(tree_file, schema='newick').length()
    
    ad = cal_gap_prop(alignment_dir, file_name) ## cal proportion of gaps in alignment somehow
    outfile = open(tree_outfile, 'a+')
    name = f.replace(trees_directory, "")
    tree_handle = open(tree_file, 'r')
    newick_string = tree_handle.readline()
    fixed_newick = fix_newick(newick_string)
    record = "".join([name,"\t",str(tl), "\t", str(ad), "\t", tree_type, "\t", fixed_newick])
    tree_handle.close()
    outfile.write(record) 
    print "Tree %s, successfully saved." %name
    outfile.close()

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



### Some important variables ###
tree_outfile = "trees.out"   # final, combined output file name
trees_directory = ""         # directory where tree files should be saved
keep_all_trees = 1           # keep all tree output files when the process is done. (0 = Delete Everything, 1 = Keep Everything)
alignment_dir = ""           # location where all alignment files are stored.
num_threads = 8              # number of threads to use in RAxML. Don't use more than physical cores in on the computer.
num_replicates = 10         # number of replicates for RAxML
tree_type = "FastTree rapid bayesian support"
###



complete_trees = done_trees(trees_directory) # get a list of the trees already built, by comparing files in working dir vs those in alignment.files
alignments = glob.glob(alignment_dir + '*.fas') # get a list of all phylip files in /alignment_dir
alignment_names = [] #names of alignment files without directory or suffix
print "\n\nGetting a list of alignment files..."
for a in alignments:
    name = a.split("/")[-1].replace(".fas", "")
    alignment_names.append(name)

s = set(complete_trees)
files_to_use = [x for x in alignment_names if x not in s] ## keep the files that are in the alignment directory, but not in the complete_trees list
print "%s files new files found.\n" % len(files_to_use)
print "-"*50+"\n"


## build the trees using RAxML
count = 1 ## status counter
for f in files_to_use:
    try:
        infile = "".join([alignment_dir, f,'.fas'])
        print "Builing tree for alignment %s, tree %s / %s..." %(f, str(count), len(files_to_use))
        subprocess.call(["FastTree", "-nt", "-gtr", "-gamma", "-quiet", "-nopr", "-out", trees_directory + f, infile])  
        write_trees(f, tree_outfile, trees_directory, tree_type)
        print "Tree %s / %s Complete." %(count, len(files_to_use))
        count += 1
    
    except:
        print "Exception for file %s, logging file for later check." %infile
        f_out = open('raxml-error.log','a')
        f_out.write(infile+"\n") # write file name to the error-log
        f_out.close() 

delete_trees(keep_all_trees, trees_directory) # delete or keep all tree files
print "-"*50
print "Tree file generation complete."
