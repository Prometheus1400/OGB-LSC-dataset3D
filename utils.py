from rdkit import Chem
import pubchempy
from numpy import loadtxt
from threading import Thread
import numpy as np
import os
import pubchempy as pcp
from tqdm import tqdm

class ThreadWithReturn(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None
    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args,
                                                **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return

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

def split_into(iterable, num_splits):
    chunk = np.array_split(np.array(iterable),num_splits)
    for i in range(len(chunk)):
        chunk[i] = list(chunk[i])
    
    return chunk

def smile2CID(smile):
    compounds = pubchempy.get_compounds(smile, namespace='smiles')
    cid = compounds[0].cid
    return cid

def CIDs2smiles(CID_list):
    Lists = []

    num_splits = len(CID_list) // 1000
    if num_splits == 0:
        num_splits = 1
    # call the pubchem api in batches of 1000
    CID_split = split_into(CID_list, num_splits)

    for cids in tqdm(CID_split):
        compounds = pcp.get_compounds(cids, as_dataframe=False)
        to_append = np.empty(len(compounds), dtype=tuple)
        for i in range(len(compounds)):
            to_append[i] = (compounds[i].canonical_smiles, compounds[i].cid)
        Lists.append(to_append)
    
    to_return = np.concatenate(Lists)
    return to_return


# def combine_databases(target_db:str, path:str) -> None:
#     """
#     combines all databases from target directory with
#     the target database
#     """
#     files = []
#     for file in os.listdir(path):
#         if file.endswith(".db") and file != target_db:
#             files.append(file)
#     writer = Writer(target_db)
#     for db in files:
#         writer.combine_with(db)