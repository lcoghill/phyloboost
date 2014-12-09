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
        
        (SPAN('Sequences'), False, URL('sequences', 'index'), [   
        (T('Browse All Sequences'), False, URL('sequences', 'index')),
        (T('Search Sequences'), False, URL('sequences', 'search'))]),
        
        (SPAN('Clusters'), False, URL('clusters', 'index'), [   
        (T('Browse All Clusters'), False, URL('clusters', 'index')),
        (T('Search Clusters'), False, URL('clusters', 'search'))]),
        
        (SPAN('Convex Subtrees'), False, URL('convexsubtrees', 'index'), [   
        (T('Browse All Convex Subtrees'), False, URL('convexsubtrees', 'index')),
        (T('Search Convex Subtrees'), False, URL('convexsubtrees', 'search'))]),

        (SPAN('Trees'), False, URL('trees', 'index'), [   
        (T('Browse All Unrooted Trees'), False, URL('trees', 'index')),
        (T('Search Unrooted Subtrees'), False, URL('trees', 'search'))]),
        
        (SPAN('About'), False, URL('about', 'index'), [   
        (T('About Phyloboost'), False, URL('about', 'index')),
        (T('Methods'), False, URL('about', 'methods')),
        (T('Wiki'), False, "https://github.com/lcoghill/phyloboost/wiki"), 
        (T('Bibliography'), False, URL('about', 'bibliography')),
        (T('Credits'), False, URL('about', 'credits'))])]
           

#if "auth" in locals(): auth.wikimenu()
