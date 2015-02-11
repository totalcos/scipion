import os
from django.conf.urls import url
import pyworkflow as pw

MEDIA_MYFIRSTMAP = os.path.join(pw.HOME, 'web', 'webtools', 'myfirstmap', 'resources')

urls = [
    (r'^resources_myfirstmap/(?P<path>.*)$', 
        'django.views.static.serve', 
        {'document_root': MEDIA_MYFIRSTMAP}
    ),
    
    url(r'^service_projects/', 'app.views_webtools.service_projects'),
    url(r'^check_project_id/$', 'app.views_webtools.check_project_id'),
    url(r'^create_service_project/$', 'app.views_webtools.create_service_project'),
    url(r'^get_testdata/$', 'app.views_webtools.get_testdata'),
    url(r'^service_content/$', 'app.views_webtools.service_content')
    
]