from django.urls import path
from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("about/", views.about, name="about"),
    path("result/", views.result, name="result"),
    path("result/donwloadcsv/", views.downloadcsv, name="downloadcsv"),
    path("result/analysis/", views.analize, name="analysis"),  
    path("result/analysis/downloadpdb/", views.downloadpdb, name="downloadpdb"),  
    path("docking/", views.docking, name="docking"),  
    path("randomforest/", views.randomforest, name="randomforest"),  
    path("shapley/", views.shapley, name="shapley"),  
    path("help/", views.help, name="help"),  
]