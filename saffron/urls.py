from django.urls import path

from . import views

urlpatterns = [
    path('', views.home, name='saffron_home'),
    path('loading/', views.loading, name='saffron_loading'),
    path('results/', views.results, name='saffron_results'),
    path('download/csv/', views.download_csv, name='saffron_download_csv'),
    path('download/excel/', views.download_excel, name='saffron_download_excel'),
    path('plot/<path:job_file>', views.plot_file, name='saffron_plot'),
    path('about/', views.about, name='saffron_about'),
]
