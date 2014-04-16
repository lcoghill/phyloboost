from StringIO import StringIO
import MySQLdb as mdb
import glob
import os

### A few important variables

input.files = "/full/path/to/blast/files/here"
model.file = "/full/path/to/model/file/here"
out.path = "/full/path/to/output/directory/here"
out.file = "file-name-here"
db.host = ""
db.user = ""
db.password = ""
db.database = ""

###


## get list of all blast result files
print "\n\nGetting a list of GI list files..."
files = glob.glob(blast.files) # get a list of all files in input.files folder
file_count = len(files)
print "%s files successfully found.\n" %file_count
print "-"*50+"\n"

## open a connection to the mysql database and set the cursor
db = mdb.connect(db.host, db.user, db.password, db.database)
cursor = db.cursor()


## get list of model organism ti values
model_list = []
model_list = [line.strip() for line in open(model.file)]

for f in files:
    original_list = []
    filtered_list = []
    print "Processing file %s..." %f
    original_list = [line.strip() for line in open(f)]
    for gi in original_list :
        sql_query = "".join(['SELECT ti FROM taxid WHERE gi=',gi,';'])
        cursor.execute(sql_query)
        ti_of_gi = cursor.fetchone()	
	if ti_of_gi :
	    ti_of_gi = str(ti_of_gi[0])
	    if ti_of_gi not in model_list:
                filtered_list.append(gi)
    out = "".join([out.path,out.file])
    out_file = open(out, 'w+')
    for item in set(filtered_list):
        print>>out_file, item
    print "Original Records: %s" %len(original_list)
    print "Unique Filtered Records: %s" %len(set(filtered_list))
    print "-"*100
db.close()
