import MySQLdb as mdb
import sys
import os
import logging
import pickle


### Read in a parse the tree file output from 'phylota-graph'

treeFile = open('tree_file.txt', 'r')

ci = [];
ci_ti = [];

for line in treeFile:
    label = line.split('	', 1 ); #split the label id from the rest of the tree
    label2 = label[0].split('_', 1 ); #split the label into ti and ci
    label2[0] = label2[0][2:]
    label2[1] = label2[1][2:]
    ci_ti.append(label2)
 

### To pull sequences from phylota mySQL db based on the ti and ci information retreived from the tree file    
config = [];
db_params = open('private/config', 'r') ## define a config file so that private db information isn't shared across code repositories. to create, renamed /private/config.example to /private/config

for line in db_params:
	line = line.rstrip('\n')
	line = line.split("=", 1)
	config.append (line[1])

database = mdb.connect(host=config[1], # your host, usually localhost
                     user=config[2], # your username
                     passwd=config[3], # your password
                     db=config[0]) # name of the data base


cluster_db = database.cursor() # define the db cursor

ti_list = []
count = 0
for cluster in ci_ti:
    #sql query to find the sequences that belong to the cluster / taxon that are present in the given tree
    #sql = "".join(["SELECT seqs.gi,seqs.seq,seqs.def FROM seqs LEFT JOIN ci_gi_184 ON seqs.gi=ci_gi_184.gi WHERE ci_gi_184.ti=2763"," AND ci_gi_184.clustid=505 AND ci_gi_184.cl_type='subtree';"])
    sql = "".join(["SELECT seqs.gi,seqs.seq,seqs.def FROM seqs LEFT JOIN ci_gi_184 ON seqs.gi=ci_gi_184.gi WHERE ci_gi_184.ti=", ci_ti[count][0]," AND ci_gi_184.clustid=",ci_ti[count][1], " AND ci_gi_184.cl_type='subtree';"])
    cluster_db.execute(sql) #execute the above sql query
    record = cluster_db.fetchall() #fetch all records that meet the query criteria and place them in the list record
    ti_list.append(record) #append the list 'record' to the tuple ti_list
    filename = "".join(["fasta/ti",ci_ti[count][0],"_ci",ci_ti[count][1],".FASTA"]) # create a string that contains the appropriate filename for each FASTA file
  

    counter = 0
    for r in ti_list:
        cur_record = str(ti_list[counter][0]) #convert the active element of the tuple ti_list to a string, and assing it to cur_record
        cur_record = cur_record.split(",") #split the cur_record string on the commas and reassign to cur_record which becomes a list
        gi = cur_record.pop(0) #pull out the gi value
        gi = gi[1:].replace("L","") #strip a leading '(' from the gi value
        seq = cur_record.pop(0).replace('\'',"")[1:]  #pull out the actual sequence, only leaving behind the elements of the description
        description = "".join(cur_record).replace('\'', "").strip(")")[1:] #join all the pieces of the description into a single string and remove problem characters and a leading space
        fasta_entry = "".join([">gi|",gi,"|",description,"\n",seq]) #parse and format strings into a single string fasta entry
        counter = counter + 1
    
    print str(count) + ": " + gi + " successful"
    count = count + 1
    
    

    