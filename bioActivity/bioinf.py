from rdkit import Chem
from rdkit.Chem import AllChem

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
    return fp_list

