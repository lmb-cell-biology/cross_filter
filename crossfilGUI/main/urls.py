from django.urls import path

from . import views

urlpatterns = [
  # path('', views.load_csv, name='load_csv'),
  path('', views.set_crossfil_page, name='set_crossfil_page'),

]

