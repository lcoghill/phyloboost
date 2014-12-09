from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from cStringIO import StringIO
from Bio import Phylo
import re


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

  rows = db().select(db.trees.ALL, limitby=(start, end))
  total_count = db.executesql('SELECT COUNT(*) FROM trees;')[0][0]
  gi_counts = {}
  for r in rows :
    gis = []
    temp_list = r['tree'].split(",")
    for t in temp_list :
      gi = re.search('\|(.+?):', t).group(1)
      if gi :
        gis.append(gi)
      gi_counts[r['id']] = len(set(gis))

  return dict(rows=rows, rec_num=rec_num, total_count=total_count, gi_counts=gi_counts)


def search():

  if request.vars.search and request.vars.field:
    rec = "".join(request.vars.search)
    field_name = "".join(request.vars.field)
    
    if field_name == "id" :
      col_name = "id"
      query = "".join(["SELECT * FROM trees WHERE id=", rec,";"])
    elif field_name == "ti_ci" :
      col_name = "ti_ci"
      q = "".join(["'%",rec,"%'"])
      query = "".join(["SELECT * FROM trees WHERE ", col_name, " LIKE ", q," limit 10;"])
    else :
      col_name = "tree"
      q = "".join(["'%gi|",rec,"%'"])
      query = "".join(["SELECT * FROM trees WHERE ", col_name, " LIKE ", q," limit 10;"])
    rows = db.executesql(query)
    
    gi_counts = {}
    for r in rows :
      gis = []
      temp_list = r[2].split(",")
      for t in temp_list :
        gi = re.search('\|(.+?):', t).group(1)
        if gi :
          gis.append(gi)
      gi_counts[r[0]] = len(set(gis))

    
    return dict(rows=rows, table_id=col_name,  gi_counts=gi_counts)
  else:
    return dict(rows=None)

def phylogram():
   
   rec = request.args(0)
   query = "".join(["SELECT * FROM trees WHERE id=", rec,";"])
   rows = db.executesql(query)
   tree = Phylo.read(StringIO(rows[0][2]), "newick")
   labels = []
 
  
   for clade in tree.find_clades() :
    if clade.name != None :
      labels.append(clade.name)

   gis = []

   for l in labels :
    gi = l.split("|")[1]
    gis.append(gi)

   return dict(rec=rec, gis=gis)

def newick():
    rec = request.args(0)
    query = "".join(["SELECT * FROM trees WHERE id=", rec,";"])
    rows = db.executesql(query)
    tree = rows[0][2]
    tree = tree.replace(",", ":1.0,")
    tree = tree.replace(")", ":1.0)")

    return dict(tree=tree)