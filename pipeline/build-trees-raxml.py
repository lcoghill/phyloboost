import subprocess
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
def write_trees(file_name, tree_outfile, trees_directory, tree_type):

    tree_file = "".join([trees_directory,"RAxML_bipartitions.", file_name])
    info_file = open("".join([trees_directory,"RAxML_info.",file_name]), "r")
    for line in info_file :
        if "Tree-Length:" in line :
            tree_length = line.split(":")[1].replace(" ", "").strip("\n")
        elif "Proportion of gaps" in line :
            ad = str(100.00 - float(line.split(": ")[1].replace("%", "").strip("\n")))
 
    info_file.close()
    outfile = open(tree_outfile, 'a+')
    print "-"*50
    print "Writing trees to file %s" %tree_outfile
    name = f.replace(trees_directory+"RAxML_bestTree.", "")
    tree_file = open(tree_file, 'r')
    newick_string = tree_file.readline()
    record = "".join([name,"\t",tree_length, "\t", ad, "\t", tree_type, "\t", newick_string])
    outfile.write(record) 
    print "Tree %s, successfully written." %name
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
num_threads = 1              # number of threads to use in RAxML. Don't use more than physical cores in on the computer.
num_replicates = 100         # number of replicates for RAxML
evo_model = "GTRGAMMA"       # model of evolution used for RAxML
p_seed = 12345               # parsimony seed for RAxML
bs_seed = 12345              # bootstrap seed for RAxML
algorithm = 'a'              # bootstrap algorithm for RAxML
tree_type = "RAxML rapid bootstrapping and subsequent ML search"
align_size_cutoff = 1000     # maximum allowed size of an alignment. ## not implemented yet.
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
        print "Builing tree for %s..." %f
        subprocess.call(["raxmlHPC-PTHREADS", "-f", "a", "-m", evo_model, "-p", str(p_seed), 
            "-x", str(bs_seed), "-#", str(num_replicates), "-T", str(num_threads), "-s", infile, "-n", f, "-w", trees_directory])
            
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
