from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC



def index():

  if request.vars.rec_count :
    rec_num = int("".join(request.vars.rec_count))

  elif request.vars.c :
    rec_num = int(request.vars.c)
  
  else:
    rec_num = 10

  if request.args(0) :
    if request.vars.d == 'next' :
      start = int(request.args(0))
      end = int(start) + rec_num

    elif request.vars.d == 'back' :
      if rec_num < int(request.args(0)) :
        start = 0
        end = int(request.args(0)) - 1
      else:
        start = (int(request.args(0)) - rec_num) - 1
        end = int(request.args(0)) - 1

    else:
      start = int(request.args(0))
      end = int(request.args(0)) + rec_num

  else:
    start = 0
    end = rec_num


  rows = db().select(db.clusters.ALL, limitby=(start, end))
  total_count = db.executesql('SELECT COUNT(*) FROM clusters;')[0][0]

  other_atts = []
  for row in rows :
    query = "".join(["SELECT gene, genome, description FROM sequences WHERE gi=", str(row['longest_gi']),";"])
    result = db.executesql(query)
    if result :
      other_atts.append(result)
    else :
      temp = [['', '', '']]
      other_atts.append(temp)

  return dict(rows=rows, rec_num=rec_num, total_count=total_count, other_atts=other_atts)

def view():

  rec = str(request.args(0))
  row = db(db.clusters.id==rec).select()
  query = "".join(["SELECT * FROM sequences WHERE gi=", str(row[0]['longest_gi']),";"])
  result = db.executesql(query)
  gis = []
  gis.append(row[0]['gi_list'].split("|"))


  return dict(row=row, result=result, gis=gis)

def fasta():
  
  rec = str(request.args(0))
  row = db(db.clusters.id==rec).select()
  records = []
  gis = row[0]['gi_list'].split("|")
  
  for gi in gis :
    seqrow = db(db.sequences.gi==gi).select().as_list()
    if seqrow :
      rec_id = "".join(["gi",str(seqrow[0]['gi']),"_ti",str(seqrow[0]['ti'])])
      record = ">" + rec_id + " " + seqrow[0]['tax_name'] + " " + seqrow[0]['description'] + "\n" + seqrow[0]['sequence'].upper()
      records.append(record)

  f_name = row[0]['pci'] + ".fas"
  
  return dict(fasta=records, f_name=f_name)

def gis():
  
  rec = str(request.args(0))
  row = db(db.clusters.id==rec).select()
  records = []
  gis = row[0]['gi_list'].split("|")
  f_name = row[0]['pci'] +".txt"
  
  return dict(fasta=gis, f_name=f_name)