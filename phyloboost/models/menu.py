# -*- coding: utf-8 -*-
# this file is released under public domain and you can use without limitations

#########################################################################
## Customize your APP title, subtitle and menus here
#########################################################################

#response.logo = A(B('web',SPAN(2),'py'),XML('&trade;&nbsp;'),
#                  _class="brand",_href="http://www.web2py.com/")
response.title = request.application.replace('_',' ').title()
response.subtitle = ''

## read more at http://dev.w3.org/html5/markup/meta.name.html
response.meta.author = 'Lyndon M Coghill <lcoghill@fieldmuseum.org>'
response.meta.keywords = ''
response.meta.generator = ''

## your http://google.com/analytics id
response.google_analytics_id = None

#########################################################################
## this is the main application menu add/remove items as required
#########################################################################
app = request.application
ctr = request.controller
response.menu = [
    (T('Home'), False, URL('default', 'index'), []),
    (T('About'), False, URL('about', 'index'), []),
    (SPAN('Sequences'), 'sequences', False, [   
        (T('Browse All Sequences'), False, URL('sequences', 'index')),
        (T('Search Sequences'), False, URL('sequences', 'search'))]),
        
        (SPAN('Clusters'), False, 'clusters', [   
        (T('Browse All Clusters'), False, URL('admin', 'default', 'site')),
        (T('Search Clusters'), False, URL('admin', 'default', 'design/%s' % app))]),
        
        (SPAN('Trees'), False, 'trees', [   
        (T('Browse All Unrooted Trees'), False, URL('trees', 'index')),
        (T('Browse All Convex Subtrees'), False, URL('convexsubtrees', 'index')),
        (T('Search Unrooted Subtrees'), False, URL('trees', 'search')),
        (T('Search Convex Subtrees'), False, URL('convexsubtrees', 'search'))])]


#if "auth" in locals(): auth.wikimenu()
