import MySQLdb as mdb
import sys
import os
import logging

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
    #ci_ti[1] = ci_ti[1][2:]
    #ci.append (ci_ti[1]) #create a list of ci values

    
print ci_ti

config = [];
db_params = open('private/config', 'r')

for line in db_params:
	line = line.rstrip('\n')
	line = line.split("=", 1)
	config.append (line[1])

#print config

database = mdb.connect(host=config[1], # your host, usually localhost
                     user=config[2], # your username
                     passwd=config[3], # your password
                     db=config[0]) # name of the data base

# you must create a Cursor object. It will let
# you execute all the query you need
cluster_db = database.cursor()


# Use all the SQL you like
ti_list = []
count = 0
for cluster in ci_ti:
    sql = "".join(["SELECT seqs.gi,seqs.seq,seqs.def FROM seqs LEFT JOIN ci_gi_184 ON seqs.gi=ci_gi_184.gi WHERE ci_gi_184.ti=", ci_ti[count][0]," AND ci_gi_184.clustid=",ci_ti[count][1], " AND ci_gi_184.cl_type='subtree';"])
    #print sql
    cluster_db.execute(sql)
    ti_list.append(cluster_db.fetchall())
    count = count + 1

print ti_list