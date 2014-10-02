# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## Customize your APP title, subtitle and menus here
#########################################################################

response.logo = A(B('web',SPAN(2),'py'),XML('&trade;&nbsp;'),
                  _class="brand",_href="http://www.web2py.com/")
response.title = request.application.replace('_',' ').title()
response.subtitle = ''

## read more at http://dev.w3.org/html5/markup/meta.name.html
response.meta.author = 'Your Name <you@example.com>'
response.meta.keywords = 'web2py, python, framework'
response.meta.generator = 'Web2py Web Framework'

## your http://google.com/analytics id
response.google_analytics_id = None

#########################################################################
## this is the main application menu add/remove items as required
#########################################################################

response.menu = [
    (T('Home'), False, URL('default', 'index'), [])
]

DEVELOPMENT_MENU = True

#########################################################################
## provide shortcuts for development. remove in production
#########################################################################

def _():
    # shortcuts
    app = request.application
    ctr = request.controller
    # useful links to internal and external resources
    response.menu += [
        (SPAN('Sequences'), 'sequences', False, [   
        (T('Browse All Sequences'), False, URL('sequences', 'view')),
        (T('Search Sequences'), False, URL('sequences', 'search'))]),
        
        (SPAN('Clusters'), False, 'clusters', [   
        (T('Browse All Clusters'), False, URL('admin', 'default', 'site')),
        (T('Search Clusters'), False, URL('admin', 'default', 'design/%s' % app))]),
        
        (SPAN('Trees'), False, 'trees', [   
        (T('Browse All Unrooted Trees'), False, URL('admin', 'default', 'site')),
        (T('Browse All Convex Subtrees'), False, URL('admin', 'default', 'design/%s' % app)),
        (T('Search Unrooted Subtrees'), False, URL('admin', 'default', 'design/%s' % app)),
        (T('Search Convex Subtrees'), False, URL('admin', 'default', 'design/%s' % app))])]
if DEVELOPMENT_MENU: _()

if "auth" in locals(): auth.wikimenu() 
