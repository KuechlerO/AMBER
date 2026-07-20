"""
URL configuration for crispr_webapp project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/6.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import path, include
from django.conf import settings


urlpatterns = [
    path('admin/', admin.site.urls),
    path('saffron/', include('saffron.urls')),
    path('', include('designer.urls')),
]

if settings.DEBUG:
    # FORCE_SCRIPT_NAME makes STATIC_URL include the subpath (/crispr-project/static/),
    # but PATH_INFO after the proxy is /static/... — so register finders at /static/.
    # Do NOT use django.conf.urls.static with a single app document_root: the first
    # /static/ mount would catch all requests and 404 files from other apps.
    urlpatterns += staticfiles_urlpatterns(prefix='/static/')
