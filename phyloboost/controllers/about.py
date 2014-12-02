def index():

    query = "".join(["SELECT * FROM versions;"])
    all_rows = db.executesql(query)
    for r in all_rows :
    	if r[3] is 1 :
    		row = r

    return dict(row=row, all_rows=all_rows)