"""crossfilGUI URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
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
from django.urls import include, path
from filebrowser.sites import site

print(site.storage.location)

urlpatterns = [
    path('admin/filebrowser/',site.urls),
    path('grappelli/', include('grappelli.urls')),
    path('', include('main.urls')),
    # path('main/', include('main.urls')), # Not included because it's redundant with main path
    path('admin/', admin.site.urls)
]

from django.core.files.storage import DefaultStorage
from filebrowser.sites import FileBrowserSite

# Default FileBrowser site
site = FileBrowserSite(name='filebrowser', storage=DefaultStorage())
