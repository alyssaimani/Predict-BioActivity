# Predict BioActivity
Predict BioActivity is a web tool designed for predicting molecular bioactivity, specifically targeting Human Acetylcholinesterase (PDB ID: 4EY7), an enzyme closely associated with Alzheimer's pathology. This web tool implements a QSAR model built using features that were extracted by generating 1024 Morgan Fingerprints representing the substructure of compounds. SHapley Additive exPlanations (SHAP) are implemented to understand locally and globally important features from the prediction results of the QSAR model. The effectiveness of the QSAR model was tested with 10-fold cross-validation, where the regression model can achieve a MAPE score of 11.10% and the classification model achieves an AUC-ROC score of 84.77%.

# How to run
 ```shell
python manage.py runserver
```
