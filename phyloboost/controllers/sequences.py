def index():

  form=FORM('Show :', INPUT(_name='count'), INPUT(_type='submit'))

  if form.vars.count :
    rec_num = form.vars.count

  else:
    rec_num = 10

  if request.args(0) :
    if request.vars.d == 'next' :
      start = int(request.args(0))
      end = int(start) + rec_num

    if request.vars.d == 'back' :
      start = (int(request.args(0)) - rec_num) - 1
      end = int(request.args(0)) - 1
    
  else:
    start = 0
    end = rec_num


  rows = db().select(db.sequences.ALL, limitby=(start, end))#.as_list()
  return dict(rows=rows, form=form, rec_num=rec_num)

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