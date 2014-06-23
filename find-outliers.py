import multiprocessing as mp
from Bio import AlignIO
import networkx as nx
import numpy as np
import glob
import itertools
import csv


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return ([e for e in t if e != None] for t in itertools.izip_longest(*args))

def csv_dump(input, f):

    split_name = f.split(".alignment")
    file_name = "".join(["csv/", split_name[0][4:], ".csv"])
    result_file = open(file_name, 'w+')
    csv_write = csv.writer(result_file, delimiter="\t")

    for item in input:
        csv_write.writerow(item)

def load_merged_nodes():
    merged_nodes = {}
    merged_file = open('ncbi/merged.dmp', 'r')
    for line in merged_file:
        line_split = line.split("\t|\t")
        new = line_split[1].replace("\t|", "").strip("\n")
        original = line_split[0]
        merged_nodes[original] = new

    return merged_nodes

def parse_seq_names(alignment):

    sequence_names = []
    ti_list = []
    for a in alignment:
        sequence_names.append(a.name)
    for n in sequence_names:
        name_split = n.split("_ti_")
        ti_list.append(name_split[1])
    return ti_list, sequence_names

def correct_merged_nodes(ti_list, merged_nodes, f):
    print "Finding any merged nodes and correcting in %s..." % f
    for ti in ti_list:
        source = str(ti)
        if source in merged_nodes:
            index = ti_list.index(ti)
            ti_list[index] = int(merged_nodes[source])
    return ti_list

def find_missing_nodes(ti_list, sequence_names, f, G):

    print "Checking that all sequence TI values in %s are in taxonomy..." % f
    missing_sequences = []
    missing_reason = []
    for ti in ti_list:
        if G.has_node(ti) is False:
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

def get_accurate_group_mean(group_dist):
	clean_dist_list = []
	within_group_cutoff = 1
	if len(clean_dist_list) == 0:	
		for d in group_dist:
			if d < within_group_cutoff + sum(int(i) for i in group_dist)/len(group_dist): # used to filter large distances before obtaining group mean.
				clean_dist_list.append(d)
		if len(clean_dist_list) == 0:	
			within_group_cutoff += 1
		else:
			return clean_dist_list

def calc_distances(ti_list, f, G):

    print "Calculating distances for %s..." % f
    mean_distances = []
    for ti in ti_list:
        group_dist = []
        for ti2 in ti_list:
            if ti != ti2:
                distance = nx.shortest_path_length(G, str(ti), str(ti2))
                group_dist.append(distance)
	clean_dist_list = get_accurate_group_mean(group_dist)
	mean_group_dist = sum(int(i) for i in clean_dist_list)/len(clean_dist_list)
        mean_distances.append(mean_group_dist)
    return mean_distances

def filter_by_distance(file_names, cutoff_threshold, merged_nodes, G):
    for f in file_names:
        print "Processing %s..." % f
        alignment = AlignIO.read(f, 'fasta')  # read in alignment file
        removed_sequences = []
        ti_list, sequence_names = parse_seq_names(alignment)  # get clean ti values and individual sequence names
        ti_list = correct_merged_nodes(ti_list, merged_nodes, f)  # correct any incorrectly labeled nodes using the merged nodes file
        bad_sequences, reason = find_missing_nodes(ti_list, sequence_names, f, G)  # are tis in taxonomy
	
	if bad_sequences:
		for seq in bad_sequences:
			bad_ti = seq.split("ti_")[1] 
			if bad_ti in ti_list:
				ti_list.remove(bad_ti)  
        
	if len(ti_list) > 1:
		mean_distances = calc_distances(ti_list, f, G)  # calc topological distances between all nodes from alignment and return
		for dist in mean_distances:
			if np.median(mean_distances)*cutoff_threshold == 0:
				cutoff = 1
           		else:
				cutoff = np.median(mean_distances)*cutoff_threshold
           	 	if dist > cutoff:
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

    # return removed_sequences
    if removed_sequences:
        csv_dump(removed_sequences, f)



'''
The cutoff_threshold is the value that will be multiplied by the total median distance for an entire alignment.
Higher values will lead to a more conservative filtering algorithm. (ie: an outlier will have be further away
taxonomically in order to be flagged.)
'''
cutoff_threshold = 3  # multiplier for mean group distance that will exclude a sequence.


graphml_file = 'ncbi.gml' # graphml file name / location
input_dir = "new-alignments/" # directory of alignment files in FASTA
numProcs = int(raw_input("Number of Cores: "))

if __name__ == '__main__':
    print "Using %s cores." %numProcs
    alignment_files = glob.glob("".join([input_dir, "*.fasta"]))  # get a list of alignment files
    print "Loading NCBI taxonomy graph..."
    G = nx.read_graphml(graphml_file) # load the graphml file into a networkx graph object
    merged_nodes = load_merged_nodes()  # get a dictionary of all nodes from merged.dmp file from ncbi
    print "Complete."
    print "-"*50
    print "\nBeginning topological distance calculations...\n"
    print "-"*50

    files_per_core = int(round(len(alignment_files) / numProcs))
    alignment_list = list(mygrouper(files_per_core, range(int(len(alignment_files)))))

    # assign the file names to the proper elements inside of the alignment_list
    # this organizes the files into a single list of lists to be processed by the mp module
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
    for c in xrange(len(alignment_list)):
        p = mp.Process(target=filter_by_distance, args=(alignment_list[c], cutoff_threshold, merged_nodes, G))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
