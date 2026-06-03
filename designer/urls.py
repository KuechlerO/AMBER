from django.urls import path
from . import views
 
urlpatterns = [
    path('', views.home, name='home'),
    path('results/', views.results, name='results'),
    path('download/excel/', views.download_excel, name='download_excel'),
    path('download/csv/', views.download_csv, name='download_csv'),
    path('loading/', views.loading, name='loading'),
    path('tutorial/', views.tutorial, name='tutorial'),
    path('about/', views.about, name='about'),
]
