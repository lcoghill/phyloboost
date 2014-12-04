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


  rows = db().select(db.sequences.ALL, limitby=(start, end))
  total_count = db.executesql('SELECT COUNT(*) FROM sequences;')[0][0]
  return dict(rows=rows, rec_num=rec_num, total_count=total_count)

def view():

  rec = str(request.args(0))
  row = db(db.sequences.id==rec).select()
  
  stats = []
  a = 0.0
  g = 0.0
  c = 0.0
  t = 0.0

  for l in row[0]['sequence']:
    if l is 'a':
      a += 1
    elif l is 'g':
      g += 1
    elif l is 'c':
      c += 1
    else:
      t += 1

  seq_len = len(row[0]['sequence'])
  stats.append(round(float(a/seq_len)*100,2))
  stats.append(round(float(g/seq_len)*100,2))
  stats.append(round(float(c/seq_len)*100,2))
  stats.append(round(float(t/seq_len)*100,2))


  return dict(row=row, stats=stats)

def search():

  #if request.vars.search and len(request.vars.search) > 3 :
  if request.vars.search and len(request.vars.search) >= 3 and request.vars.field:
    rec = "".join(request.vars.search)
    q = "".join(["'",rec,"%'"])
    col_name = "".join(request.vars.field)
    #rows = db(db.sequences.tax_name.like(q)).select(limitby=(0, 5))
    query = "".join(["SELECT * FROM sequences WHERE ", col_name, " LIKE ", q," limit 10;"])
    #SELECT * FROM sequences WHERE tax_name LIKE 'Sus%' limit 10;
    rows = db.executesql(query)
    return dict(rows=rows, table_id=col_name)
  else:
    return dict(rows=None)

def fasta():
  rec = request.args(0)
  row = db(db.sequences.id==rec).select().as_list()
  rec_id = "".join(["gi",str(row[0]['gi']),"_ti",str(row[0]['ti'])])
  record = SeqRecord(Seq(row[0]['sequence'].upper(),
                   IUPAC.protein),
                   id=rec_id, name=row[0]['tax_name'],
                   description=row[0]['description'])

  f_name = rec_id+".fas"
  
  return dict(fasta=unicode(record.format("fasta"), "utf-8"), f_name=f_name)