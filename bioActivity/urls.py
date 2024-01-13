from django.urls import path
from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("about/", views.about, name="about"),
    path("result/", views.result, name="result"),
    path("result/donwload/", views.download, name="download"),
    path("result/analysis/", views.analize, name="analysis"),  
    path("docking/", views.docking, name="docking"),  
    path("randomforest/", views.randomforest, name="randomforest"),  
    path("shapley/", views.shapley, name="shapley"),  
]