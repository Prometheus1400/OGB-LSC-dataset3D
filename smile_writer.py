from collections import namedtuple
from database_api import Reader, Writer
import numpy as np
from utils import split_into, ThreadWithReturn
import time
from random import shuffle
from tqdm import tqdm, trange
import pandas as pd
import pubchempy as pcp
from rdkit import Chem

def CIDs2smiles(CID_list):
    Lists = []

    num_splits = len(CID_list) // 10
    if num_splits == 0:
        num_splits = 1
    # call the pubchem api in batches of 1000
    CID_split = split_into(CID_list, num_splits)

    for cids in (CID_split):
        compounds = pcp.get_compounds(cids, as_dataframe=False)
        to_append = np.empty(len(compounds), dtype=tuple)
        for i in range(len(compounds)):
            to_append[i] = (compounds[i].canonical_smiles, compounds[i].cid)
        if len(compounds) != 0:
            Lists.append(to_append)
    if len(Lists) != 0:
        return np.concatenate(Lists)
    return Lists

def smiles2CIDs(SMILES_list):
    Lists = []

    num_splits = len(SMILES_list) // 10
    if num_splits == 0:
        num_splits = 1
    # call the pubchem api in batches of 1000
    smile_split = split_into(SMILES_list, num_splits)

    for smiles in smile_split:
        print(smiles)
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        # compounds = pcp.get_properties(['cid'], smile, 'smiles', as_dataframe=True)
        to_append = np.empty(len(compounds), dtype=tuple)
        for i in range(len(compounds)):
            to_append[i] = (compounds[i].canonical_smiles, compounds[i].cid)
        if len(compounds) != 0:
            Lists.append(to_append)
    if len(Lists) != 0:
        return np.concatenate(Lists)
    return Lists

reader = Reader("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
smiles = reader.get_all_SMILES()
del(reader)

new_smiles = np.empty(len(smiles), dtype=tuple)
writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
for i in trange(len(smiles)):
    try:
        new_smile = Chem.MolToSmiles(Chem.MolFromSmiles(smiles[i]))
    except Exception:
        new_smile = smiles[i]
        print(f"failed to convert {smiles[i]}")
    new_smiles[i] = (new_smile, smiles[i])

writer.add_smiles(new_smiles)

# smiles = list(pd.read_csv("/data3/kaleb.dickerson2001/Datasets/OGB-LSC-3D/pcqm4m_kddcup2021/raw/data.csv.gz")['smiles'].to_numpy())[:10]
# print(smiles)
# print(smiles2CIDs(smiles))


# if __name__ == "__main__":
#     NUM_THREADS = 2
#     total_split = split_into(new_total, NUM_THREADS)

#     T = []
#     for i in range(NUM_THREADS):
#         T.append(ThreadWithReturn(target=CIDs2smiles, args=(total_split[i],)))
#         T[i].start()
    
#     arr = np.empty(NUM_THREADS, dtype=object)
#     for i,t in enumerate(T):
#         arr[i] = t.join()

#     smiles_with_CIDS = np.concatenate(arr)
    
#     time.sleep(30)
#     writer = Writer("/data3/kaleb.dickerson2001/Datasets/PubChem3D/coord_db.db")
#     writer.add_smiles(smiles_with_CIDS, type='canonical')