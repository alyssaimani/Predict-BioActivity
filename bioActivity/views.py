from django.shortcuts import render, HttpResponse
from .forms import InsertFile
import pandas as pd


# Create your views here.
def home(request):
    if request.method == "POST":
        
        form = InsertFile(request.POST, request.FILES)
        
        if form.is_valid(): 
            file = form.cleaned_data['file']
            predict(file)  
 
    else:
        form = InsertFile()
    return render(request, 'home.html', {'form': form})

def about(request):
    return render(request, 'about.html')

def predict(file_path):
    df = pd.read_csv(file_path, delimiter=',')
    print(df.values)