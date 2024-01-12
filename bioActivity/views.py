from django.shortcuts import render, HttpResponse, redirect
from .forms import InsertFile
from .utils import get_fingerprints, get_graph
from .models import SmilesData, CurrSmilesData
import pandas as pd
from PIL import Image 
import pickle, os, csv, shap
import dill as pickle_dill
import matplotlib.pyplot as plt
from django.http import HttpResponse
from rdkit.Chem import Draw
from rdkit import Chem
from django.shortcuts import render

def predict(file):
    ext = os.path.splitext(file.name)[1]
    if ext.lower() == '.csv':
        df_smiles = pd.read_csv(file, delimiter=',', header=None)
        df_smiles = [smiles[0] for smiles in df_smiles.values]
    elif ext.lower() == '.xlsx':
        df_smiles = pd.read_excel(file, header=None)
        df_smiles = [smiles[0] for smiles in df_smiles.values]
    elif ext.lower() == '.txt':
        file.seek(0)
        df_smiles = file.read().decode('utf-8')
        df_smiles = df_smiles.splitlines()

    reg_model = pickle.load(open("assets/models/regression_model.pkl", 'rb'))
    class_model = pickle.load(open("assets/models/classification_model.pkl", 'rb'))
    fingerprints,_ = get_fingerprints(df_smiles)
    pIC50_pred = reg_model.predict(fingerprints)
    class_pred = class_model.predict(fingerprints)
    run_SHAP(fingerprints, class_model)
    results = pd.DataFrame()
    results['smiles'] = df_smiles
    results['pIC50'] =  pIC50_pred
    results['bio_class'] = [True if pred == 1 else False for pred in class_pred]
    return results  

def run_SHAP(fingeprints, class_model):
    # SHAP visualization 
    df_fingerprints = pd.DataFrame(fingeprints)
    explainer = shap.TreeExplainer(class_model)
    shap_values = explainer.shap_values(df_fingerprints, check_additivity=False)
    plt.figure()
    shap.summary_plot(shap_values, df_fingerprints, feature_names=df_fingerprints.columns, show=False)
    # Save the plot to a file
    plot_path = 'static/images/shap_summary.png'
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()
    return

def save_data(data):
    CurrSmilesData.objects.all().delete()
    for _,row in data.iterrows():
        SmilesData.objects.get_or_create(smiles=row['smiles'], pic50=row['pIC50'], bio_class=row['bio_class'])
        CurrSmilesData.objects.create(smiles=row['smiles'], pic50=row['pIC50'], bio_class=row['bio_class'])

    return

def download(request):
    smiles_data = SmilesData.objects.all()
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="smiles_data.csv"'
    writer = csv.writer(response)
    writer.writerow(['smiles', 'pIC50', 'bioactivity_class'])
    for data in smiles_data:
        writer.writerow([data.smiles, data.pic50, data.bio_class])

    return response

# Create your views here.
def home(request):

    if request.method == "POST":  
        form = InsertFile(request.POST, request.FILES)
        if form.is_valid(): 
            file = form.cleaned_data['file']
            results = predict(file)  
            results = results.drop_duplicates(subset=['smiles'])
            save_data(results)
            results = results.to_dict(orient='records')
            return redirect('result')
    else:
        with open('assets/pdb_files/clean_receptor.pdb', 'r') as file:
            protein_pdb = file.read()
        form = InsertFile()

    return render(request, 'home.html', {'form':form, 'protein_pdb': protein_pdb})

def result(request):
    results = CurrSmilesData.objects.all()
    return render(request, 'result.html', {'results':results})

def about(request):
    return render(request, 'about.html')

def analize(request):
    smiles = request.GET.get('smiles', '')
    pic50 = request.GET.get('pIC50', '')
    bio_class = request.GET.get('bio_class', '')
    fingerprints,bit_list = get_fingerprints([smiles])
    mfpvector = [index for index, value in enumerate(fingerprints[0]) if value != 0]
    fp_images = []
    molecule = Chem.MolFromSmiles(smiles)
    for idx in mfpvector:
        svg_image = Draw.DrawMorganBit(molecule, idx, bit_list[0])
        fp_images.append(svg_image) 
    images_vectors = zip(fp_images, mfpvector)
    
    # Lime visualization
    class_model = pickle.load(open("assets/models/classification_model.pkl", 'rb'))
    explainer = pickle_dill.load(open('assets/models/lime_explainer.pkl', 'rb'))
    df_fingerprints = pd.DataFrame(fingerprints)
    instance = df_fingerprints.iloc[0]
    explanation = explainer.explain_instance(instance, class_model.predict_proba)
    
    # Convert to HTML
    lime_html = explanation.as_html()   
    return render(request, 'analysis.html', {'smiles':smiles, 'pic50':pic50, 'bio_class':bio_class, 'images_vectors':images_vectors, 'lime_html':lime_html})

def docking(request):
    with open('assets/pdb_files/clean_receptor.pdb', 'r') as file:
        receptor = file.read()

    with open('assets/pdb_files/highest.pdb', 'r') as file:
        highest = file.read()

    with open('assets/pdb_files/lowest.pdb', 'r') as file:
        lowest = file.read()
    
    return render(request, 'docking.html', {'receptor':receptor, 'highest':highest})
