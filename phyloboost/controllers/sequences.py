def index():
  theme = "smoothness"
  for x in (
      # 'DataTables-1.8.1/media/js/jquery.js',  # there's already a newer (Bootstrap-compatible) jQuery loaded!
      'DataTables-1.8.1/media/js/jquery.dataTables.min.js',
      'DataTables-1.8.1/media/css/bootstrap_table.css',
      'DataTables-1.8.1/media/ui/css/%s/jquery-ui-1.8.5.custom.css' % theme):
      response.files.append(URL('static',x))

  colnames = ["Id", "TI", "Taxon", "Cluster ID",
              "Sequence", "Genome", "Gene", "Accession",
              "GI", "Version", "Mol. Type", "GB Div", "Description"]
  widths = ["2%", "2%", "5%", "2%", "10%", "5%", "5%",
            "5%", "2%", "2%", "5%", "2%", "53%"]
  tid = "sequences"
  table = TABLE(_id=tid, _class="display")
  table.append(THEAD(TR(*[ TH(f, _width=w)
                           for f, w in zip(colnames, widths) ])))
  table.append(TBODY(TR(TD("Loading data from server",
                           _colspan=len(colnames),
                           _class="dataTables_empty"))))
  table.append(TFOOT(TR(
      TH(INPUT(_name="search_id",
               _style="width:100%",_class="search_init",
               _title="search Id" )),
      TH(INPUT(_name="search_focal_clade",
               _style="width:100%",_class="search_init",
               _title="search focal clade" )),
      TH(INPUT(_name="search_study",
               _style="width:100%",_class="search_init",
               _title="search study" )),
      TH(INPUT(_name="search_type",
               _style="width:100%",_class="search_init",
               _title="search tree type" )),
      TH(),
      TH(INPUT(_name="search_uploaded",
               _style="width:100%",_class="search_init",
               _title="search uploaded" )),
      TH(INPUT(_name="search_person",
               _style="width:100%",_class="search_init",
               _title="search person" )))))

  
  rows = db(db.sequences.id > 10).select()
  return dict(tid=tid, rows=rows)
