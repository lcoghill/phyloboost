def index():

    query = "".join(["SELECT * FROM versions;"])
    all_rows = db.executesql(query)

    if request.args(0) :
        for r in all_rows :
            if r[1] == str(request.args(0)) :
                row = r

    else :
        row = all_rows[-1]
        
    return dict(row=row, all_rows=all_rows)

def methods():

    query = "".join(["SELECT * FROM versions;"])
    all_rows = db.executesql(query)
    for r in all_rows :
        if r[3] is 1 :
            row = r

    return dict(row=row, all_rows=all_rows)

def bibliography():

    query = "".join(["SELECT * FROM versions;"])
    all_rows = db.executesql(query)
    for r in all_rows :
        if r[3] is 1 :
            row = r

    return dict(row=row, all_rows=all_rows)

def credits():

    query = "".join(["SELECT * FROM versions;"])
    all_rows = db.executesql(query)
    for r in all_rows :
        if r[3] is 1 :
            row = r

    return dict(row=row, all_rows=all_rows)