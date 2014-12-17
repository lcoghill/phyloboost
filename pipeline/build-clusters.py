import MySQLdb as mdb



input_file = "phylota_184_trees.tre"
out_dir = "clusters/"
db_host = "localhost"
db_user = "root"
db_password = "reelab14"
db_database = "phylota"



### Read in a parse the tree file output from trees file

treefile = open(input_file, 'r')
raw_id = []
ids = []
for t in treefile:
    tree_id = t.split("\t")
    tree_id = tree_id[0].split("_")
    tree_id[0] = tree_id[0][2:]
    tree_id[1] = tree_id[1][2:]
    ids.append(tree_id)

## connect to the database and set a cursor
database = mdb.connect(host = db_host,  # your host, usually localhost
                     user = db_user,  # your username
                     passwd = db_password,  # your password
                     db = db_database)  # name of the data base
cluster_db = database.cursor()  # define the db cursor


count = 0  # status count
for i in ids:
    # sql query to find the sequences that belong to the cluster / taxon that are present in the given tree
    sql = "".join(["SELECT seqs.gi,seqs.seq FROM seqs LEFT JOIN ci_gi_184 ON seqs.gi=ci_gi_184.gi WHERE ci_gi_184.ti=", ids[count][0], " AND ci_gi_184.clustid = ", ids[count][1], " AND ci_gi_184.cl_type='subtree';"])
    cluster_db.execute(sql)  # execute the above sql query
    record = cluster_db.fetchall()  # fetch all records that meet the query criteria and place them in the list record
    record = list(record)  # convert tuple to list
    # create a string that contains the appropriate filename for each FASTA file
    filename = "".join([out_dir, "ti", ids[count][0], "_ci", ids[count][1], ".fas"])
    f = open(filename, 'w+')
    
    for r in record:
        cur_record = list(r)  # convert this element of the record list from tuple to a list
        cur_record[0] = str(cur_record[0])  # convert the GID value from long to string
        # join all elements of the cur_record list into a single string, added formatting for FASTA style
        cur_record = "".join([">gi|", str(cur_record[0]), "\n", cur_record[1]])
        f.write("%s\n" % cur_record)    
    
    print "Sequences for tree " + str(count) + " successfully written to " + filename + "."
    count += 1
    f.close()    
