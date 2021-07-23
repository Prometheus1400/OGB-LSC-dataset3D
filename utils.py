from rdkit import Chem
import pubchempy
from numpy import loadtxt


def CID_from_SMILES(smile):
    compounds = pubchempy.get_compounds(smile, namespace='smiles')
    cid = compounds[0].cid
    smiles_1 = compounds[0].canonical_smiles

    mol = Chem.MolFromSmiles(smile)
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True)


    mol_1 = Chem.MolFromSmiles(smiles_1)
    smiles_1_canonical = Chem.MolToSmiles(mol_1, canonical=True)

    if smiles_1_canonical == smiles_canonical:
        # print('cid: ', cid)
        # print(smiles_canonical)
        return cid
    else:
        # print('Not found')
        # print(smiles_1_canonical)
        # print(smiles_canonical)
        return -1

def get_coordinates(file_path, error_log=None) -> str:
    """
    input: .inp file containing optimized coords
    return: string of coords seperated by space
    """
    try:
        data = loadtxt(file_path,delimiter="\n",dtype=str)[9:][:-1]
        length = len(data)
        # to_return = empty(length, dtype=tuple)
        to_return = ""
        first = False
        for i in range(length):
            line = data[i].split()
            s = " ".join(line[2:])
            to_return += " " + s
        return to_return[1:] #remove leading space
    except:
        if error_log != None:
            error_log.write(f"error in get_coordinates at file: {file_path}\n")
        return -1

# def get_smiles(file_path) -> str:
#     """
#     input: .info file
#     return: smiles string
#     """
#     # simply returns smiles
#     try:
#         smile = loadtxt(file_path, dtype=str)[3]
#         return smile
#     except:
#         error_log.write(f"error in get_smiles at file: {file_path}\n")
#         return "NONE"