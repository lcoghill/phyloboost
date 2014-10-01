from Bio import SeqIO
import MySQLdb
import glob
import re


def fetch_gi(acc, conn, cur, acc_gi_table):

    cur.execute("SELECT gi FROM %s WHERE acc=%s" %(acc_gi_table, acc))
    row = cur.fetchone()
    if row is not None:
        gi = row[0]
    else:
        gi = '00000'

    return str(gi)

def fetch_taxid(gi, conn, cur, gi_taxid_table):

    cur.execute("SELECT taxid FROM %s WHERE gi=%s" %(gi_taxid_table, gi))
    row = cur.fetchone()
    
    if row is not None:
        taxid = row[0]
    else: 
        taxid = "00000"
        

    return str(taxid)

def convert_date(date):

    date_vars = date.split("-")
    monthDict={'JAN':'01', 'FEB':'02', 'MAR':'03', 'APR':'04', 'MAY':'05', 'JUN':'06', 'JUL':'07', 'AUG':'08', 'SEP':'09', 'OCT':'10', 'NOV':'11', 'DEC':'12'}
    day = date_vars[0]
    month = monthDict[date_vars[1]]
    year = date_vars[2]
    new_date = year + "-" + month + "-" + day

    return new_date
    


## variables that need set
db_user = ''
db_passwd = ''
db = ''
genbank_dir = ''
gi_taxid_table = ''
acc_gi_table = ''



conn = MySQLdb.connect(user=db_user, passwd=db_passwd, db=db)
cur = conn.cursor()
genbank_file_list = glob.glob(genbank_dir + '*.seq')
csv_f = open('genbank_eukrayotes.csv', 'a')


for gbfile in genbank_file_list:
    handle = open(gbfile, 'rU')
    print "Converting Genbank file %s..." %gbfile
    for record in SeqIO.parse(handle, "gb") :
        acc = record.id
        description = record.description.replace(",", "").strip(".")
        seq = record.seq
        gb_div = re.search('gb(.+?).seq', gbfile).group(1).upper()
    

        if record.features[0].qualifiers['organism']:
            organism =  record.features[0].qualifiers['organism'][0]
        elif record.annotations['organism']:
            organism = record.annotations['organism'][0]
        else:
            organism = 'ignotus'


        if record.features[0].qualifiers['mol_type']:
            mol_type = "".join(record.features[0].qualifiers['mol_type'])
        else:
            mol_type = 'ignotus'


        if record.annotations['taxonomy']:
            taxonomy = "|".join(record.annotations['taxonomy'])
        else:
            taxonomy = 'ignotus'


        if record.annotations['date']:
            date = convert_date(record.annotations['date'])
        else:
            date = '00-00-0000'


        if record.annotations['gi']:
            gi = record.annotations['gi']
        else:
            gi = fetch_gi(acc, conn, cur)

    
        line = ",".join([gi, mol_type, organism, taxonomy, date, gb_div, description, acc, str(seq)])
        csv_f.write(line+'\n')
    handle.close()

conn.close()
csv_f.close()
