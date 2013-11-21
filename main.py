import MySQLdb
import os
import logging
import library.genbank
logging.basicConfig(level=logging.INFO)
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

### A few import static variables
library.genbank.email = 'me@my.address.com'


db = MySQLdb.connect(host="localhost", # your host, usually localhost
                     user="root", # your username
                      passwd="reelab", # your password
                      db="phylota") # name of the data base

# you must create a Cursor object. It will let
#  you execute all the query you need
ti = db.cursor() 


# Use all the SQL you like
ti.execute("SELECT * FROM ci_gi_184 WHERE ti=54375 and clustid=19")

gi_list = []
# print all the first cell of all the rows
for row in ti.fetchall():
    gi_list.append(str(row[3]))
    

print gi_list
genes = library.genbank.fetch_genelist(gi_list)
print genes