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
    path('results/screen-plot/status/', views.screen_plot_status, name='screen_plot_status'),
    path('results/screen-plot/warmup/', views.screen_plot_warmup, name='screen_plot_warmup'),
    path('results/screen-plot/overview/', views.screen_plot_overview, name='screen_plot_overview'),
    path('results/screen-plot/prewarm-full/', views.screen_plot_prewarm_full, name='screen_plot_prewarm_full'),
    path('results/screen-plot/full-status/', views.screen_plot_full_status, name='screen_plot_full_status'),
    path('results/screen-plot/', views.screen_enrichment_plot, name='screen_enrichment_plot'),
    path('results/structure/<str:accession>.pdb', views.structure_pdb, name='structure_pdb'),
]
