from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC



def index():

  if request.vars.rec_count :
    rec_num = int(request.vars.rec_count)
    start = int(request.vars.row1) -1
    end = start + rec_num

  else :  
    if request.vars.r :
      rec_num = int(request.vars.r)
  
    else:
      rec_num = 10

    if request.args(0) :
      if request.vars.d == 'next' :
        start = int(request.args(0))
        end = start + rec_num

      elif request.vars.d == 'back' :
        if int(request.args(0)) < rec_num :
          start = 0
          end = rec_num
        else :
          end = int(request.args(0)) - 1
          start = end - rec_num
      else :
        start = int(request.args(0))
        end = start + rec_num
        
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

def search():

  if request.vars.search and request.vars.field:
    rec = "".join(request.vars.search)
    field_name = "".join(request.vars.field)
    
    if field_name == "id" :
      col_name = "id"
      query = "".join(["SELECT * FROM clusters WHERE id=", rec,";"])
    elif field_name == "gi" :
      col_name = "gi_list"
      q = "".join(["'%",rec,"%'"])
      query = "".join(["SELECT * FROM clusters WHERE ", col_name, " LIKE ", q," limit 10;"])
    else :
      col_name = "pci"
      q = "".join(["'%ti",rec,"%'"])
      query = "".join(["SELECT * FROM clusters WHERE ", col_name, " LIKE ", q," limit 10;"])
    rows = db.executesql(query)
    
    other_atts = []
    for row in rows :
      query = "".join(["SELECT gene, genome, description FROM sequences WHERE gi=", str(row[4]),";"])
      result = db.executesql(query)
      if result :
        other_atts.append(result)
      else :
        temp = [['', '', '']]
        other_atts.append(temp)
    
    return dict(rows=rows, table_id=col_name, other_atts=other_atts)
  else:
    return dict(rows=None)