from Bio import AlignIO
from datetime import datetime
from ete2 import Tree
import numpy as np
import glob
import multiprocessing as mp
import itertools
import csv

#  functions


def csv_dump(input, f):

    split_name = f.split(".alignment")
    file_name = "".join(["csv/", split_name[0][4:], ".csv"])
    result_file = open(file_name, 'w+')
    csv_write = csv.writer(result_file, delimiter="\t")

    for item in input:
        csv_write.writerow(item)

def parse_seq_names(alignment):

    sequence_names = []
    ti_list = []
    for a in alignment:
        sequence_names.append(a.name)
    for n in sequence_names:
        name_split = n.split("_ti_")
        ti_list.append(name_split[1])
    return ti_list, sequence_names


def get_all_nodes(tree):

    print "Getting all nodes from the tree..."
    all_tis = []
    for node in tree.traverse():
        all_tis.append(int(node.name))
    print "Complete."

    return all_tis


def load_merged_nodes():
    merged_nodes = {}
    merged_file = open('ncbi/merged.dmp', 'r')
    for line in merged_file:
        line_split = line.split("\t|\t")
        new = line_split[1].replace("\t|", "").strip("\n")
        original = line_split[0]
        merged_nodes[original] = new

    return merged_nodes


def correct_merged_nodes(ti_list, merged_nodes, f):
    print "Finding any merged nodes and correcting in %s..." % f
    for ti in ti_list:
        source = str(ti)
        if source in merged_nodes:
            index = ti_list.index(ti)
            ti_list[index] = int(merged_nodes[source])
    return ti_list


def find_missing_nodes(ti_list, sequence_names, all_tis, f):

    print "Checking that all sequence TI values in %s are in taxonomy..." % f
    missing_sequences = []
    missing_reason = []
    for ti in ti_list:
        if int(ti) not in all_tis:
            bad_index = ti_list.index(ti)
            missing_sequences.append(sequence_names[bad_index])
            missing_reason.append("TI not found in taxonomy")

    del_nodes_file = open('ncbi/delnodes.dmp', 'r')
    del_nodes = []
    for line in del_nodes_file:
        del_nodes.append(int(line.replace("\t|", "").strip("\n")))

    for ti in ti_list:
        if int(ti) in del_nodes:
            bad_index = ti_list.index(ti)
            bad_seq = sequence_names[bad_index]
            if bad_seq not in missing_sequences:
                missing_sequences.append(sequence_names[bad_index])
                missing_reason.append("TI found in deleted nodes file")

    return missing_sequences, missing_reason


def calc_distances(ti_list, f, taxo_tree):

    print "Calculating distances for %s..." % f
    mean_distances = []
    for ti in ti_list:
        group_dist = []
        for ti2 in ti_list:
            if ti != ti2:
                distance = taxo_tree.get_distance(str(ti), target2=str(ti2), topology_only=True)
                group_dist.append(distance)

        mean_group_dist = np.mean(group_dist)
        mean_distances.append(int(mean_group_dist))

    return mean_distances


def filter_by_distance(file_names, cutoff_threshold, merged_nodes, all_tis, taxo_tree):
    for f in file_names:
        print "Processing %s..." % f
        alignment = AlignIO.read(f, 'fasta')  # read in alignment file
        removed_sequences = []
        ti_list, sequence_names = parse_seq_names(alignment)  # get clean ti values and individual sequence names
        ti_list = correct_merged_nodes(ti_list, merged_nodes, f)  # correct any incorrectly labeled nodes using the merged nodes file
        bad_sequences, reason = find_missing_nodes(ti_list, sequence_names, all_tis, f)  # are tis in taxonomy
        mean_distances = calc_distances(ti_list, f, taxo_tree)  # calc topological distances between all nodes from alignment and return
        for dist in mean_distances:
            cutoff = np.median(mean_distances)*cutoff_threshold
            if dist >= cutoff:
                bad_index = mean_distances.index(dist)
                seq_name = sequence_names[bad_index]
                if seq_name not in bad_sequences:
                    bad_sequences.append(sequence_names[bad_index])
                    reason_string = "".join(["Exceeds Cutoff: ", str(cutoff), " with Distance: ", str(dist)])
                    reason.append(reason_string)

        for seq in bad_sequences:
            temp_list = []
            index = bad_sequences.index(seq)
            temp_list.append(f[5:])
            temp_list.append(seq)
            temp_list.append(reason[index])
            removed_sequences.append(temp_list)

        #return removed_sequences
        if removed_sequences:
            csv_dump(removed_sequences, f)

def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return ([e for e in t if e != None] for t in itertools.izip_longest(*args))

# end of function definitions


startTime = datetime.now()
cutoff_threshold = 2  # multiplier for mean group distance that will exclude a sequence.
numProcs = int(raw_input("Number of Cores: "))


if __name__ == '__main__':
    print "Using %s cores." %numProcs

    # input alignment file parsing
    input_dir = "test/"
    alignment_files = glob.glob("".join([input_dir, "*.alignment.fasta"]))  # get a list of alignment files
    print "Loading NCBI taxonomy..."
    taxo_tree = Tree('ncbi_taxonomy.newick', format=1)  # import ncbi taxonomy newick string
    print "complete."
    all_tis = get_all_nodes(taxo_tree)  # get list of all nodes in ncbi taxonomy tree
    merged_nodes = load_merged_nodes()  # get a dictionary of all nodes from merged.dmp file from ncbi
    print "-"*50
    print "\nBeginning topological distance calculations...\n"
    print "-"*50

    files_per_core = int(round(len(alignment_files) / numProcs))
    alignment_list = list(mygrouper(files_per_core, range(int(len(alignment_files)))))

# assign the file names to the proper elements inside of the alignment_list
    track_alignments = 0
    core_count = 0
    for l in alignment_list:
        count = 0
        while count < len(alignment_list[core_count]):
            alignment_list[core_count][count] = alignment_files[track_alignments]
            track_alignments += 1
            count += 1

        core_count += 1


# Split job into appropriate number of cores
    processes = []
    cycle_count = 1
    for c in xrange(len(alignment_list)):
        print "Processing file %s / %s..." % (cycle_count, len(alignment_files))
        cycle_count += 1
        p = mp.Process(target=filter_by_distance, args=(alignment_list[c], cutoff_threshold, merged_nodes, all_tis, taxo_tree))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

print "Script execution took %s" % (datetime.now()-startTime)