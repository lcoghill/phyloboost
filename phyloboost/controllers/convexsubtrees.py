from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from cStringIO import StringIO
from Bio import Phylo


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

  rows = db().select(db.convex_subtrees.ALL, limitby=(start, end))
  total_count = db.executesql('SELECT COUNT(*) FROM convex_subtrees;')[0][0]

  return dict(rows=rows, rec_num=rec_num, total_count=total_count)



def view():

  rec = str(request.args(0))
  row = db(db.convex_subtrees.id==rec).select()
  tree = Phylo.read(StringIO(row[0]['tree']), "newick")
  labels = []
 
  for clade in tree.find_clades() :
    if clade.name != None :
      labels.append(clade.name)

  gis = []
  tis = []

  for l in labels :
    gis.append(l.split("_")[0][2:])
    tis.append(l.split("_")[1][2:])

  gis = set(gis)
  tis = set(tis)
  return dict(row=row, gis=gis, tis=tis)

def search():

  if request.vars.search and request.vars.field:
    rec = "".join(request.vars.search)
    field_name = "".join(request.vars.field)
    
    if field_name == "id" :
      col_name = "id"
      query = "".join(["SELECT * FROM convex_subtrees WHERE id=", rec,";"])
    elif field_name == "root_taxon" :
      col_name = "root_taxon"
      q = "".join(["'",rec,"%'"])
      query = "".join(["SELECT * FROM convex_subtrees WHERE ", col_name, " LIKE ", q," limit 10;"])
    elif field_name == "root_ti" :
      col_name = "root_ti"
      q = "".join(["'",rec,"%'"])
      query = "".join(["SELECT * FROM convex_subtrees WHERE ", col_name, " LIKE ", q," limit 10;"])
    elif field_name == "gi" :
      col_name = "tree"
      q = "".join(["'%gi",rec,"%'"])
      query = "".join(["SELECT * FROM convex_subtrees WHERE ", col_name, " LIKE ", q," limit 10;"])
    else :
      col_name = "tree"
      q = "".join(["'%ti",rec,"%'"])
      query = "".join(["SELECT * FROM convex_subtrees WHERE ", col_name, " LIKE ", q," limit 10;"])
      
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
  
  return dict(fasta=record.format("fasta"), f_name=f_name)

def phylogram():
   
   rec = str(request.args(0))
   row = db(db.convex_subtrees.id==rec).select()
   tree = Phylo.read(StringIO(row[0]['tree']), "newick")
   labels = []
 
  
   for clade in tree.find_clades() :
    if clade.name != None :
      labels.append(clade.name)

   gis = []
   tis = []

   for l in labels :
    gis.append(l.split("_")[0][2:])
    tis.append(l.split("_")[1][2:])

   gis = set(gis)
   tis = set(tis)

   return dict(rec=rec, row=row, gis=gis, tis=tis)

def newick():
    rec = request.args(0)
    query = "".join(["SELECT * FROM convex_subtrees WHERE id=", rec,";"])
    rows = db.executesql(query)
    tree = rows[0][4]
    tree = rows[0][4]
    tree = tree.replace(",", ":1.0,")
    tree = tree.replace(")", ":1.0)")

    return dict(tree=tree)