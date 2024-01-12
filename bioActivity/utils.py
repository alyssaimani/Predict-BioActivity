from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import base64
from PIL import Image
from io import BytesIO

def get_fingerprints(smiles_list):
    radius = 2
    nBits = 1024
    bit_list = []
    morgan_fp = []
    molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    for mol in molecules:
        bit = {}
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,useChirality=True, radius=radius, nBits = nBits, bitInfo=bit)
        morgan_fp.append(fingerprint)
        bit_list.append(bit)
    fp_list = [list(fingerprint) for fingerprint in morgan_fp]
    return fp_list, bit_list

def get_graph():
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64decode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    return graph

