import MySQLdb as mdb
import sys
import os
import logging
import pickle
import itertools

def parse_config():
    config_params = [];
    config = open('private/config', 'r') 
    for line in (line for line in config if not line.startswith('###')):
        line = line.rstrip('\n')
        line = line.split("=")
        config_params.append (line[1])

    return config_params

### Read in a parse the tree file output from trees file

treefile = open('trees.txt', 'r')
raw_id = []
id = []
for t in treefile:
    tree_id = t.split("\t")
    tree_id = tree_id[0].split("_")
    tree_id[0] = tree_id[0][2:]
    tree_id[1] = tree_id[1][2:]
    id.append(tree_id)

### To pull sequences from phylota mySQL db based on the ti and ci information retreived from the tree file    
## Gather import configuration information from the config file
params = parse_config() # retreive the params from the config file


database = mdb.connect(host=params[1], # your host, usually localhost
                     user=params[2], # your username
                     passwd=params[3], # your password
                     db=params[0]) # name of the data base
cluster_db = database.cursor() # define the db cursor

count = 0
for i in id:
    #sql query to find the sequences that belong to the cluster / taxon that are present in the given tree
    sql = "".join(["SELECT seqs.gi,seqs.seq FROM seqs LEFT JOIN ci_gi_184 ON seqs.gi=ci_gi_184.gi WHERE ci_gi_184.ti=", id[count][0]," AND ci_gi_184.clustid=",id[count][1], " AND ci_gi_184.cl_type='subtree';"])
    cluster_db.execute(sql) #execute the above sql query
    record = cluster_db.fetchall() #fetch all records that meet the query criteria and place them in the list record
    record = list(record) #convert tuple to list
    filename = "".join(["output/fasta/ti",id[count][0],"_ci",id[count][1],".fasta"]) # create a string that contains the appropriate filename for each FASTA file
    f = open(filename,'w+')
    
    for r in record:
        cur_record = list(r) #convert this element of the record list from tuple to a list
        cur_record[0] = str(cur_record[0]) #convert the GID value from long to string
        cur_record = "".join([">gi|",str(cur_record[0]),"\n",cur_record[1]]) #join all elements of the cur_record list into a single string, added formatting for FASTA style
        f.write("%s\n" % cur_record)    
    
    print "Sequences for tree " + str(count) + " successfully written to " + filename + "."
    count = count + 1
    f.close()    
