from django.shortcuts import render, HttpResponse
from .forms import InsertFile
from .bioinf import get_fingerprints
import pandas as pd
import pickle



# Create your views here.
def home(request):
    if request.method == "POST":
        
        form = InsertFile(request.POST, request.FILES)
        
        if form.is_valid(): 
            file = form.cleaned_data['file']
            results = predict(file)  
            results = results.to_dict(orient='records')
            return render(request, 'result.html', {'results':results})
    else:
        form = InsertFile()
    return render(request, 'home.html', {'form':form})

def about(request):
    return render(request, 'about.html')

def predict(file_path):
    df_smiles = pd.read_csv(file_path, delimiter=',')
    reg_model = pickle.load(open("assets/regression_model.pkl", 'rb'))
    class_model = pickle.load(open("assets/classification_model.pkl", 'rb'))
    df_smiles = [smiles[0] for smiles in df_smiles.values]
    fingerprints = get_fingerprints(df_smiles)
    pIC50_pred = reg_model.predict(fingerprints)
    class_pred = class_model.predict(fingerprints)
    results = pd.DataFrame()
    results['smiles'] = df_smiles
    results['pIC50'] = pIC50_pred
    results['class'] = class_pred
    return results

def download():
    return