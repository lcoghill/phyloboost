import datetime, socket, os, sys
import MySQLdb
from gluon.dal import MySQLAdapter
from ConfigParser import SafeConfigParser
from gluon.tools import *

MySQLAdapter.driver = MySQLdb
defaults = dict(host="localhost", user="guest", password="guest", dbname="phyloboost")
conf = SafeConfigParser(defaults)
user = password = dbname = host = ''


if os.path.isfile("applications/%s/private/localconfig" % request.application):
    conf.read("applications/%s/private/localconfig" % request.application)
    host = conf.get("db", "host")
    user = conf.get("db", "user")
    password = conf.get("db", "password")
    dbname = conf.get("db", "dbname")

else:
    conf.read("applications/%s/private/config" % request.application)
    host = conf.get("db", "host")
    user = conf.get("db", "user")
    password = conf.get("db", "password")
    dbname = conf.get("db", "dbname")

db = DAL("mysql://%s:%s@%s/%s" % (user, password, host, dbname), migrate=False, migrate_enabled=False )

## mail = Mail()                    # mailer
auth = Auth(globals(),db)           # authentication/authorization
crud = Crud(globals(),db)           # for CRUD helpers using auth
service = Service(globals())        # for json, xml, jsonrpc, xmlrpc, amfrpc
plugins = PluginManager()
host = socket.gethostname()
crud.settings.auth = None           # =auth to enforce authorization on crud
auth.define_tables( migrate=False )                           # creates all needed tables

try:
    host = conf.get("hosting", "hostdomain")
except:
    host = "localhost:8000"

db.define_table('convex_subtrees', 
    Field('id', 'id'),                       
    Field('ci', 'integer'),                  
    Field('date', 'date'),                   
    Field('tree', 'text'),                   
    Field('method', 'string'),               
    )

db.define_table('sequences', 
    Field('id', 'id'),                       
    Field('ti', 'integer'),                  
    Field('tax_name', 'string'),             
    Field('phylota_ci', 'integer'),          
    Field('sequence', 'text'),               
    Field('genome', 'string'),               
    Field('gene', 'text'),                   
    Field('accession', 'string'),            
    Field('gi', 'integer'),                  
    Field('pb_version', 'integer'),          
    Field('mol_type', 'text'),               
    Field('genbank_div', 'string'),          
    Field('description', 'text'),            
    )

db.define_table('trees', 
    Field('id', 'id'),                       
    Field('ci', 'integer'),                  
    Field('ti', 'integer'),                  
    Field('tree', 'text'),                   
    Field('method', 'string'),               
    Field('date', 'date'),                   
    )

db.define_table('versions', 
    Field('id', 'id'),                       
    Field('version', 'integer'),             
    Field('creation_date', 'date'),          
    Field('ncbi_taxonomy', 'date'),          
    Field('phylota_version', 'integer'),     
    Field('genbank_version', 'integer'),     
    Field('repbase_version', 'integer'),     
    )
